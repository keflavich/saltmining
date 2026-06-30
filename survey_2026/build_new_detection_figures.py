"""Build per-source diagnostic moment-map panels for the 'first reported
here' sources in Table 4.

For each source we render a single multi-panel PNG combining the brightest
NaCl, KCl, H2O 232 GHz, and SiS/SO mom0 maps (up to 4 panels). The PNGs
are saved next to the source's analysis_products dir under a `figures/`
subdir for inclusion in main.tex.

Targets are determined from data_summary's footnote system (any row with
* footnote, i.e., new this work).
"""
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ImageNormalize, AsymmetricPercentileInterval

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
OUT_FIG = Path("/orange/adamginsburg/salt/demography_2026/figures")
OUT_FIG.mkdir(parents=True, exist_ok=True)

# Map display label (Table-4 row) -> directory name on disk.
NEW_DETECTION_TARGETS = [
    ("IRAS 15412-5359",  "G326.6618+00.5207"),
    ("IRAS 18174-1612",  "G015.0357-00.6795"),
    ("I17016-4124",      "G345.5043+00.3480"),
    ("IRAS 17233-3606",  "IRAS17233-3606"),
]

# Per source, choose up to 4 mom0 maps to render. Priority order.
LINE_PRIORITIES = [
    ("NaCl", r"NaCl"),
    ("KCl",  r"KCl"),
    ("H2O 232", r"H2O.*232"),
    ("SiS",  r"SiS"),
    ("SO",   r"SO_"),
    ("H2O", r"H2O_"),
]


import pandas as pd  # noqa: E402

_TARGET_COORDS = None


def _target_coord(name):
    global _TARGET_COORDS
    if _TARGET_COORDS is None:
        df = pd.read_csv(Path("/orange/adamginsburg/salt/survey_2026/data/sources_L4_d2.csv"))
        _TARGET_COORDS = {r["name"]: (float(r["ra_deg"]), float(r["dec_deg"]))
                          for _, r in df.iterrows()}
    # Aliases for hyphen/underscore variants
    return (_TARGET_COORDS.get(name)
            or _TARGET_COORDS.get(name.replace("_", "-"))
            or _TARGET_COORDS.get(name.replace("-", "_")))


def _sep_arcsec(ra1, dec1, ra2, dec2):
    dra = (ra1 - ra2) * np.cos(np.radians(dec2)) * 3600.0
    ddec = (dec1 - dec2) * 3600.0
    return float(np.sqrt(dra * dra + ddec * ddec))


def best_proposal_dir(target_dir):
    """Pick the proposal subdir whose brightest on-target source is closest
    to the catalog target coord (within 5"), tie-break by line-measurement
    count. Falls back to plain max-lines if NO proposal has an on-target
    source (i.e., the target has no matched analysis)."""
    cat = _target_coord(target_dir.name)
    candidates = []
    for p in sorted(target_dir.glob("2*")):
        meas = p / "line_measurements.csv"
        if not meas.exists():
            continue
        try:
            n_lines = sum(1 for _ in meas.open()) - 1
        except OSError:
            n_lines = 0
        bid_info = brightest_on_target_source(p, cat)
        if bid_info is None:
            continue
        bid, sep = bid_info
        candidates.append((sep, -n_lines, p, bid))
    if not candidates:
        return None
    candidates.sort()  # smallest sep first, then most lines
    sep0, _, pdir, _ = candidates[0]
    if sep0 > 5.0:
        return None
    return pdir


def brightest_on_target_source(prop_dir, cat_coord):
    """Return (src_id, sep_arcsec) for the brightest continuum source within
    5\" of cat_coord. Falls back to the absolute brightest source if no
    on-target match exists."""
    cont = prop_dir / "continuum_sources.csv"
    if not cont.exists():
        return None
    try:
        df = pd.read_csv(cont)
    except pd.errors.EmptyDataError:
        return None
    if df.empty or "peak_Jybeam" not in df.columns:
        return None
    if cat_coord is None:
        bid = int(df.loc[df["peak_Jybeam"].idxmax(), "id"])
        return bid, np.inf
    df = df.copy()
    df["sep"] = df.apply(lambda r: _sep_arcsec(float(r["ra_deg"]),
                                                 float(r["dec_deg"]),
                                                 *cat_coord), axis=1)
    on = df[df["sep"] <= 5.0]
    if on.empty:
        return None
    row = on.sort_values("peak_Jybeam", ascending=False).iloc[0]
    return int(row["id"]), float(row["sep"])


def brightest_source(prop_dir):
    cat_coord = _target_coord(prop_dir.parent.name)
    info = brightest_on_target_source(prop_dir, cat_coord)
    return info[0] if info else None


def find_mom0(prop_dir, src_id, line_label, pattern):
    """Return the highest-priority mom0 fits matching pattern at src_id."""
    import re
    src_dir = prop_dir / f"source_{src_id:02d}"
    if not src_dir.is_dir():
        return None
    rx = re.compile(pattern, re.IGNORECASE)
    candidates = [p for p in src_dir.glob("*_mom0.fits") if rx.search(p.name)]
    # Sort: prefer v=0 (ground state, less risk of confusion) for NaCl/KCl,
    # else first.
    def key(p):
        n = p.name
        v0 = 0 if "_v0_" in n else 1
        return (v0, n)
    candidates.sort(key=key)
    return candidates[0] if candidates else None


def render_panel(ax, fits_path, title):
    with fits.open(fits_path) as hdul:
        data = hdul[0].data
        hdr = hdul[0].header
    if data is None:
        ax.set_axis_off()
        ax.set_title(f"{title} (no data)")
        return
    while data.ndim > 2:
        data = data[0]
    wcs = WCS(hdr).celestial
    norm = ImageNormalize(data, interval=AsymmetricPercentileInterval(2, 99.5))
    im = ax.imshow(data, origin="lower", cmap="inferno", norm=norm)
    ax.set_title(title, fontsize=10)
    # Strip ticks for compactness
    ax.tick_params(axis="both", which="both", labelsize=7)


def build_figure(display, target):
    tdir = ANALYSIS / target
    if not tdir.is_dir():
        for alias in (target.replace("-", "_"), target.replace("_", "-")):
            if (ANALYSIS / alias).is_dir():
                tdir = ANALYSIS / alias
                break
    if not tdir.is_dir():
        print(f"  no analysis dir for {target}")
        return None
    pd_ = best_proposal_dir(tdir)
    if pd_ is None:
        print(f"  no proposal dir under {tdir}")
        return None
    bid = brightest_source(pd_)
    if bid is None:
        print(f"  no brightest source in {pd_}")
        return None
    panels = []
    for label, pat in LINE_PRIORITIES:
        m = find_mom0(pd_, bid, label, pat)
        if m is not None:
            panels.append((label, m))
        if len(panels) == 4:
            break
    if not panels:
        print(f"  no mom0 maps for {target} src_{bid:02d}")
        return None
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4), squeeze=False)
    for ax, (label, fp) in zip(axes[0], panels):
        render_panel(ax, fp, f"{label}: {fp.stem.replace('source_'+f'{bid:02d}_','')}")
    fig.suptitle(f"{display} (src_{bid:02d}, {pd_.name})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    safe = target.replace("/", "_").replace(" ", "_")
    outp = OUT_FIG / f"newdet_{safe}_mom0.png"
    fig.savefig(outp, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {outp.name}  ({n} panels)")
    return outp


def main():
    for display, target in NEW_DETECTION_TARGETS:
        print(f"== {display}  ({target})")
        build_figure(display, target)


if __name__ == "__main__":
    main()
