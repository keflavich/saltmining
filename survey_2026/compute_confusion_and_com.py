"""Compute two diagnostic columns for data_summary.tex:

f(>5sigma)
----------
For each (target, proposal): scan all cube files used by line_pipeline (from
line_measurements.csv "cube" column), restrict to channels whose rest frequency
is in the 219.2--220.8 GHz canonical hot-core window, and compute the fraction
of those channels (averaged over the brightest source's small aperture) that
exceed 5 * sigma_native. High value -> line confusion / hot-core spectrum.

COM detection flag
------------------
True if any complex-organic line (CH3OH, CH3CN, HC3N, CH3OCHO, C2H5CN,
CH3OCH3, NH2CHO, t-HCOOH) was detected at >=5 sigma at the brightest source
in line_measurements.csv.

Writes:
  data/data_summary_aux.csv  with columns target, proposal, f5sigma, com
"""
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy import wcs

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
UVDIR = ROOT / "uvdata"
OUT = ROOT / "data/data_summary_aux.csv"

SIGMA_THRESH = 5.0
COM_FRAC_THRESHOLD = 0.05  # f(>5sigma) above this => COM-rich/hot-core

COM_PATTERNS = [
    re.compile(p, re.IGNORECASE) for p in [
        r"^CH3OH",      # methanol
        r"^CH3CN",      # methyl cyanide
        r"^HC3N",       # cyanoacetylene
        r"^CH3OCHO",    # methyl formate
        r"^C2H5CN",     # ethyl cyanide
        r"^CH3OCH3",    # dimethyl ether
        r"^NH2CHO",     # formamide
        r"^t-HCOOH",    # trans-formic
        r"^HCOOH",
        r"^HNCO",       # isocyanic
    ]
]


def brightest_source_id(cont_csv: Path):
    if not cont_csv.exists() or cont_csv.stat().st_size == 0:
        return None
    try:
        df = pd.read_csv(cont_csv)
    except pd.errors.EmptyDataError:
        return None
    if df.empty or "peak_Jybeam" not in df.columns:
        return None
    return int(df.loc[df["peak_Jybeam"].idxmax(), "id"])


def src_radec(cont_csv: Path, bid: int):
    """Return (ra_deg, dec_deg) for the source. Cube-grid pixel coordinates
    differ from continuum-image grid; the position must be re-projected via
    WCS per-cube rather than reused from continuum_sources.csv x/y."""
    df = pd.read_csv(cont_csv)
    row = df[df["id"] == bid].iloc[0]
    return float(row["ra_deg"]), float(row["dec_deg"])


def fraction_above_5sigma(cube_path: Path, ra_deg: float,
                            dec_deg: float) -> float | None:
    """Read a small spectrum around the WCS-projected (ra, dec) on this cube,
    compute fraction of channels above 5 * mad_std."""
    if not cube_path.exists():
        return None
    try:
        with fits.open(cube_path, memmap=True) as hdul:
            hdr = hdul[0].header
            data = hdul[0].data
    except (OSError, fits.VerifyError):
        return None
    if data is None or data.ndim < 3:
        return None
    arr = np.squeeze(data)
    if arr.ndim != 3:
        return None
    nz, ny, nx = arr.shape
    try:
        w = wcs.WCS(hdr).celestial
    except (wcs.InvalidTransformError, ValueError):
        return None
    # Project source RA/Dec onto this cube's celestial WCS
    try:
        xp, yp = w.wcs_world2pix(ra_deg, dec_deg, 0)
        x = int(round(float(xp)))
        y = int(round(float(yp)))
    except (ValueError, TypeError):
        return None
    if not (0 <= x < nx and 0 <= y < ny):
        return None
    y0, y1 = max(0, y - 1), min(ny, y + 2)
    x0, x1 = max(0, x - 1), min(nx, x + 2)
    spec = np.nanmean(arr[:, y0:y1, x0:x1], axis=(1, 2))
    if not np.isfinite(spec).any():
        return None
    from astropy.stats import mad_std
    sigma = mad_std(spec, ignore_nan=True)
    if not np.isfinite(sigma) or sigma <= 0:
        return None
    n_above = int(np.sum(spec > SIGMA_THRESH * sigma))
    return float(n_above) / int(np.isfinite(spec).sum())


def is_com_line(line: str) -> bool:
    return any(p.match(line) for p in COM_PATTERNS)


def main():
    import sys
    targets_filter = set(sys.argv[1:]) if len(sys.argv) > 1 else None
    rows = []
    # Preserve existing rows for targets we're not recomputing.
    existing = {}
    if targets_filter and OUT.exists():
        try:
            df_prev = pd.read_csv(OUT)
            for _, r in df_prev.iterrows():
                key = (str(r["target"]), str(r["proposal"]))
                existing[key] = r.to_dict()
        except pd.errors.EmptyDataError:
            pass
    for proposal_dir in sorted(ANALYSIS.glob("*/2*")):
        target = proposal_dir.parent.name
        proposal = proposal_dir.name
        if targets_filter and target not in targets_filter:
            # Keep pre-existing row unchanged
            if (target, proposal) in existing:
                rows.append(existing[(target, proposal)])
            continue
        meas = proposal_dir / "line_measurements.csv"
        if not meas.exists():
            continue
        try:
            df = pd.read_csv(meas)
        except pd.errors.EmptyDataError:
            continue
        if df.empty or "snr" not in df.columns:
            continue
        bid = brightest_source_id(proposal_dir / "continuum_sources.csv")
        if bid is None:
            continue
        sub = df[df["source"] == bid]
        # f(>5sigma) over the whole observed band
        ra_deg, dec_deg = src_radec(proposal_dir / "continuum_sources.csv", bid)
        fracs = []
        for cube_name in sub["cube"].dropna().unique():
            cube_path = UVDIR / proposal / target / cube_name
            f5 = fraction_above_5sigma(cube_path, ra_deg, dec_deg)
            if f5 is not None:
                fracs.append(f5)
        f5sigma = float(np.mean(fracs)) if fracs else None
        # COM proxy: explicit COM line detected, OR confusion fraction above
        # COM_FRAC_THRESHOLD (line-crowding => hot-core spectrum)
        com_lines = [L for L in sub["line"].astype(str) if is_com_line(L)]
        com_explicit = bool(sub[(sub["line"].isin(com_lines)) & (sub["snr"] >= 5)].shape[0])
        com_det = com_explicit or (f5sigma is not None and f5sigma >= COM_FRAC_THRESHOLD)
        rows.append(dict(target=target, proposal=proposal,
                          f5sigma=f5sigma, com=com_det))
        print(f"[{target} / {proposal}] f={f5sigma!r}  com={com_det}  cubes={len(fracs)}",
              flush=True)

    df_out = pd.DataFrame(rows)
    df_out.to_csv(OUT, index=False)
    print(f"\nwrote {OUT} ({len(df_out)} rows)")


if __name__ == "__main__":
    main()
