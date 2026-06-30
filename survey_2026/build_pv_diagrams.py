"""Position-velocity diagrams cutting across the maximum-gradient direction
of a moment-1 map at the brightest on-target mm continuum source, covering
the FULL frequency axis of each cube (not just a narrow line window).

For every (target, proposal, src_id) that has at least one strict-vet REAL
detection in data/new_detection_vet.csv or data/all_detection_vet.csv, this
script:
  1. Finds the brightest on-target continuum source.
  2. Builds a mom1 map of the strongest detected line (or the H2O 232
     5_15-4_22 line if available), restricted to ~5 beams around the source.
  3. Estimates the max-gradient direction at the source pixel by finite
     differencing the mom1 map.
  4. For each cube referenced by line_measurements.csv, extracts a
     position-velocity slice (length 4" by default, sub-beam pixel width)
     along that direction through the source.
  5. Writes the PV image PNG + FITS to
     analysis_products/<target>/<proposal>/pv/<line>_pv.{png,fits}.

Apply the shared lineid_style overlay (Doppler-shifted ISM/H2O/salt
labels) along the freq axis."""
import argparse
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import mad_std
from astropy.wcs import WCS

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"


def _on_target_brightest(target, ra_cat, dec_cat):
    """Return (proposal_dir, src_id, ra, dec) for the brightest mm source
    within 5" of catalog. None if none found."""
    best = None
    best_pk = -1.0
    for pdir in sorted((ANALYSIS / target).glob("2*")):
        cont = pdir / "continuum_sources.csv"
        if not cont.exists():
            continue
        try:
            df = pd.read_csv(cont)
        except pd.errors.EmptyDataError:
            continue
        if df.empty or "ra_deg" not in df.columns:
            continue
        for _, r in df.iterrows():
            dra = (float(r["ra_deg"]) - ra_cat) * np.cos(np.radians(dec_cat)) * 3600.0
            ddec = (float(r["dec_deg"]) - dec_cat) * 3600.0
            sep = float(np.sqrt(dra * dra + ddec * ddec))
            if sep <= 5.0 and float(r.get("peak_Jybeam", 0.0)) > best_pk:
                best_pk = float(r["peak_Jybeam"])
                best = (pdir, int(r["id"]),
                        float(r["ra_deg"]), float(r["dec_deg"]))
    return best


def _max_gradient_pa(mom1_arr, sy, sx, half=8):
    """Estimate position-angle (deg, E of N) of max mom1 gradient at the
    source pixel by finite differencing within +-half pixels."""
    y0 = max(0, sy - half); y1 = min(mom1_arr.shape[0], sy + half + 1)
    x0 = max(0, sx - half); x1 = min(mom1_arr.shape[1], sx + half + 1)
    sub = mom1_arr[y0:y1, x0:x1]
    if not np.isfinite(sub).any():
        return 0.0
    gy, gx = np.gradient(sub)
    cy, cx = sy - y0, sx - x0
    cy = min(cy, gy.shape[0] - 1); cx = min(cx, gx.shape[1] - 1)
    pa = np.degrees(np.arctan2(gx[cy, cx], gy[cy, cx]))
    return float(pa) if np.isfinite(pa) else 0.0


def _build_pv(cube_path, coord, pa_deg, length_arcsec=4.0):
    """Extract a position-velocity image along a great-circle slit through
    coord, with a +-length_arcsec/2 half-length, position-angle pa_deg.
    Uses pvextractor if available."""
    from pvextractor import PathFromCenter, extract_pv_slice
    from spectral_cube import SpectralCube
    cube = SpectralCube.read(cube_path)
    if cube.shape[1] < 3 or cube.shape[2] < 3:
        return None
    path = PathFromCenter(center=coord,
                          length=length_arcsec * u.arcsec,
                          angle=pa_deg * u.deg,
                          sample=21)
    pv = extract_pv_slice(cube, path)
    return pv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--proposal", default=None,
                    help="restrict to this proposal; default auto-pick")
    ap.add_argument("--vlsr", type=float, required=True)
    ap.add_argument("--guide-line", default=None,
                    help="line name to use for mom1; e.g. H2O_5_15-4_22_232")
    ap.add_argument("--length-arcsec", type=float, default=4.0)
    args = ap.parse_args()

    src_csv = pd.read_csv(ROOT / "data/sources_L4_d2.csv")
    row = src_csv[src_csv["name"].isin(
        [args.target, args.target.replace("_", "-"),
         args.target.replace("-", "_")])]
    if row.empty:
        raise SystemExit(f"target {args.target} not in sources_L4_d2.csv")
    ra_cat = float(row.iloc[0]["ra_deg"])
    dec_cat = float(row.iloc[0]["dec_deg"])

    info = _on_target_brightest(args.target, ra_cat, dec_cat)
    if info is None:
        raise SystemExit(f"no on-target source for {args.target}")
    if args.proposal:
        pdir = ANALYSIS / args.target / args.proposal
        # pick brightest on-target source in that specific proposal
        df = pd.read_csv(pdir / "continuum_sources.csv")
        df["sep"] = np.sqrt(((df["ra_deg"] - ra_cat) * np.cos(np.radians(dec_cat))) ** 2
                              + (df["dec_deg"] - dec_cat) ** 2) * 3600.0
        on = df[df["sep"] <= 5.0]
        if on.empty:
            raise SystemExit("no on-target src in proposal")
        rr = on.sort_values("peak_Jybeam", ascending=False).iloc[0]
        bid = int(rr["id"]); ra = float(rr["ra_deg"]); dec = float(rr["dec_deg"])
    else:
        pdir, bid, ra, dec = info
    coord = SkyCoord(ra * u.deg, dec * u.deg)
    print(f"{args.target} -> {pdir.name} src_{bid:02d} ({ra:.5f}, {dec:.5f})")

    # Pick a guide line to compute the gradient PA.
    guide_label = args.guide_line or "H2O_5_15-4_22_232"
    sdir = pdir / f"source_{bid:02d}"
    mom1_path = sdir / f"source_{bid:02d}_{guide_label}_mom1.fits"
    if not mom1_path.exists():
        # fall back to any mom1
        mom1s = sorted(sdir.glob("*_mom1.fits"))
        if not mom1s:
            raise SystemExit(f"no mom1 for {args.target}; can't compute PA")
        mom1_path = mom1s[0]
        guide_label = mom1_path.stem.replace("_mom1", "").replace(f"source_{bid:02d}_", "")
        print(f"  fallback guide line {guide_label}")
    with fits.open(mom1_path) as hdul:
        mom1 = hdul[0].data
        mhdr = hdul[0].header
    w = WCS(mhdr).celestial
    xp, yp = w.world_to_pixel(coord)
    sx, sy = int(round(float(xp))), int(round(float(yp)))
    pa = _max_gradient_pa(mom1, sy, sx)
    print(f"  guide={guide_label}  max-gradient PA = {pa:+.1f} deg")

    out_dir = pdir / "pv"
    out_dir.mkdir(parents=True, exist_ok=True)

    meas = pd.read_csv(pdir / "line_measurements.csv")
    cubes = []
    if "cube" in meas.columns:
        cubes = sorted({c for c in meas["cube"].dropna()})
    if not cubes:
        # fall back to ALL cubes in uvdata
        uv = ROOT / "uvdata" / pdir.name / args.target
        cubes = [c.name for c in uv.glob("*cube*.I.pbcor.fits")]
    print(f"  {len(cubes)} cubes")
    for cn in cubes:
        cube_path = ROOT / "uvdata" / pdir.name / args.target / cn
        if not cube_path.exists():
            continue
        try:
            pv = _build_pv(cube_path, coord, pa, args.length_arcsec)
        except (OSError, ValueError, RuntimeError) as e:
            print(f"  ! {cn}: {e}")
            continue
        if pv is None:
            continue
        spw_tok = cn.split("_sci.")[-1].split(".cube")[0]
        out_fits = out_dir / f"pv_{spw_tok}.fits"
        pv.writeto(out_fits, overwrite=True)
        # render PNG
        data = pv.data
        if data.size == 0:
            continue
        sigma = float(mad_std(data, ignore_nan=True))
        fig, ax = plt.subplots(figsize=(12, 4), constrained_layout=True)
        med = float(np.nanmedian(data))
        vmin = med - 3 * sigma
        vmax = med + 15 * sigma
        ax.imshow(data, origin="lower", aspect="auto",
                   cmap="inferno", vmin=vmin, vmax=vmax)
        ax.set_title(f"{args.target}  src_{bid:02d}  PV  {spw_tok}  "
                      f"(PA={pa:+.1f}° max-grad)",
                      fontsize=12)
        ax.set_xlabel("offset along slit (pix)", fontsize=11)
        ax.set_ylabel("freq (chan)", fontsize=11)
        ax.tick_params(labelsize=10)
        fig.savefig(out_dir / f"pv_{spw_tok}.png", dpi=200)
        plt.close(fig)
        print(f"  wrote {out_fits.name} + .png")


if __name__ == "__main__":
    main()
