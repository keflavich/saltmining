"""Cut down ALMA pbcor cubes to small boxes around candidate sources.

For each uvdata/<proposal>/<target>/<member>*.cube*.I.pbcor.fits cube:
  - Look up continuum_sources.csv for that (target, proposal).
  - Build a cutout bounding box = the union of per-source boxes, where
    each per-source box has half-width = max(2000 AU at the source
    distance, 3*beam_FWHM).
  - Write the cutout to uvdata/<proposal>/<target>/cutouts/<orig>.cut.fits
    (preserves full spectral axis, drops empty stokes).

Idempotent: skips if the cutout already exists with at least the same
NAXIS1/2 as the requested box.

Usage:
  python cube_cutout_2000au.py                # dry-run; just report sizes
  python cube_cutout_2000au.py --apply        # actually write cutouts
  python cube_cutout_2000au.py --apply --target G326.6618+00.5207
                                              # restrict to one target
"""
import argparse
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import astropy.units as u
import astropy.coordinates as ac
from astropy.io import fits
from astropy.wcs import WCS

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
UVDIR = ROOT / "uvdata"
ANALYSIS = ROOT / "analysis_products"
SRC_CSV = ROOT / "data/sources_L4_d2.csv"

CUBE_RE = re.compile(r"\.cube.*\.I\.pbcor\.fits$", re.IGNORECASE)
TARGET_RADIUS_AU = 2000.0
MIN_BEAM_FACTOR = 3.0


def srcs_for_field(target, proposal):
    p = ANALYSIS / target / proposal / "continuum_sources.csv"
    if not p.exists() or p.stat().st_size == 0:
        return None
    try:
        return pd.read_csv(p)
    except pd.errors.EmptyDataError:
        return None


def cube_files(target_dir):
    return [p for p in target_dir.iterdir()
             if CUBE_RE.search(p.name)]


def coord_to_pix(wcs, ra, dec):
    sky = ac.SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    try:
        x, y = wcs.world_to_pixel(sky)
    except (ValueError, TypeError):
        return None
    return float(x), float(y)


def cube_radius_pix(h, dist_kpc):
    """Cutout half-width in pixels: max(2000 AU, 3 beam_FWHM) at this cube's
    pixel scale."""
    pix_arcsec = abs(float(h.get("CDELT1", 1.0))) * 3600.0
    if pix_arcsec <= 0:
        return None
    if dist_kpc and dist_kpc > 0:
        au_arcsec = TARGET_RADIUS_AU / (dist_kpc * 1000.0)
    else:
        au_arcsec = 3.0
    bmaj = h.get("BMAJ") or 0.0
    beam_arcsec = float(bmaj) * 3600.0
    half_arcsec = max(au_arcsec, MIN_BEAM_FACTOR * beam_arcsec)
    return max(2, int(np.ceil(half_arcsec / pix_arcsec)))


def cutout_box(h, srcs, dist_kpc, target_ra=None, target_dec=None):
    """Return (x0, y0, x1, y1) bounding-box in pixels, or None if the box
    spans (almost) the whole cube — meaning no point in cutting."""
    wcs = WCS(h).celestial
    nx, ny = int(h["NAXIS1"]), int(h["NAXIS2"])
    half = cube_radius_pix(h, dist_kpc)
    if half is None:
        return None
    boxes = []
    # Per-source boxes (from continuum_sources.csv)
    if srcs is not None and not srcs.empty:
        for _, r in srcs.iterrows():
            px = coord_to_pix(wcs, float(r["ra_deg"]), float(r["dec_deg"]))
            if px is None:
                continue
            x, y = px
            if not (0 <= x < nx and 0 <= y < ny):
                continue
            boxes.append((x - half, y - half, x + half, y + half))
    # Fallback: target coordinate (may differ from any continuum source)
    if target_ra is not None and target_dec is not None:
        px = coord_to_pix(wcs, target_ra, target_dec)
        if px is not None:
            x, y = px
            if 0 <= x < nx and 0 <= y < ny:
                boxes.append((x - half, y - half, x + half, y + half))
    if not boxes:
        return None
    x0 = max(0, int(np.floor(min(b[0] for b in boxes))))
    y0 = max(0, int(np.floor(min(b[1] for b in boxes))))
    x1 = min(nx, int(np.ceil(max(b[2] for b in boxes))))
    y1 = min(ny, int(np.ceil(max(b[3] for b in boxes))))
    if x1 - x0 < 4 or y1 - y0 < 4:
        return None
    if (x1 - x0) * (y1 - y0) >= 0.85 * nx * ny:
        # Less than 15% reduction — not worth cutting
        return None
    return x0, y0, x1, y1


def update_wcs(h, x0, y0, x1, y1):
    """Patch CRPIX1/CRPIX2 and NAXIS1/NAXIS2 for a (x0:x1, y0:y1) cutout."""
    h2 = h.copy()
    h2["NAXIS1"] = int(x1 - x0)
    h2["NAXIS2"] = int(y1 - y0)
    if "CRPIX1" in h2:
        h2["CRPIX1"] = float(h2["CRPIX1"]) - x0
    if "CRPIX2" in h2:
        h2["CRPIX2"] = float(h2["CRPIX2"]) - y0
    return h2


def maybe_cut(cube_path, target, proposal, srcs, dist_kpc, target_ra,
              target_dec, apply=False):
    """Returns (orig_bytes, cut_bytes, status)."""
    with fits.open(cube_path, memmap=True) as hdul:
        h = hdul[0].header.copy()
        ndim = int(h.get("NAXIS", 0))
        nx = int(h.get("NAXIS1", 0))
        ny = int(h.get("NAXIS2", 0))
        nz = int(h.get("NAXIS3", 0))
        if ndim < 3 or nx == 0 or ny == 0:
            return (cube_path.stat().st_size, None, "skip-2d")
    box = cutout_box(h, srcs, dist_kpc, target_ra, target_dec)
    if box is None:
        return (cube_path.stat().st_size, None, "skip-cover")
    x0, y0, x1, y1 = box
    outdir = cube_path.parent / "cutouts"
    out = outdir / (cube_path.stem + ".cut.fits")
    orig_bytes = cube_path.stat().st_size
    nx_new, ny_new = x1 - x0, y1 - y0
    est_bytes = int(orig_bytes * (nx_new * ny_new) / (nx * ny))
    if not apply:
        return (orig_bytes, est_bytes, f"dryrun-{nx_new}x{ny_new}")
    if out.exists() and out.stat().st_size >= 0.9 * est_bytes:
        return (orig_bytes, out.stat().st_size, "already")
    outdir.mkdir(exist_ok=True)
    # Stream-write: open original, slice, write out
    with fits.open(cube_path, memmap=True) as hdul:
        data = hdul[0].data  # (chan, y, x) or (stokes, chan, y, x)
        if data.ndim == 4:
            cut = data[:, :, y0:y1, x0:x1]
        elif data.ndim == 3:
            cut = data[:, y0:y1, x0:x1]
        else:
            return (orig_bytes, None, "skip-shape")
        new_h = update_wcs(h, x0, y0, x1, y1)
        # Drop history about old cube size
        for k in ("DATAMAX", "DATAMIN"):
            if k in new_h:
                del new_h[k]
        fits.PrimaryHDU(data=cut, header=new_h).writeto(
            out, overwrite=True, output_verify="fix"
        )
        # Preserve auxiliary HDUs (e.g., BEAMS table); rewrite using fits.append
        for hdu in hdul[1:]:
            with fits.open(out, mode="append", memmap=True) as out_hdul:
                out_hdul.append(hdu.copy())
    return (orig_bytes, out.stat().st_size, "wrote")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true",
                    help="actually write cutouts (default: dry-run)")
    ap.add_argument("--target", default=None,
                    help="restrict to this target")
    ap.add_argument("--proposal", default=None,
                    help="restrict to this proposal")
    args = ap.parse_args()

    src_df = pd.read_csv(SRC_CSV).set_index("name")
    total_in = 0
    total_out = 0
    n_files = 0
    n_wrote = 0
    n_skip = 0
    for prop_dir in sorted(UVDIR.iterdir()):
        if not prop_dir.is_dir():
            continue
        if args.proposal and prop_dir.name != args.proposal:
            continue
        for target_dir in sorted(prop_dir.iterdir()):
            if not target_dir.is_dir():
                continue
            target = target_dir.name
            if args.target and target != args.target:
                continue
            cubes = cube_files(target_dir)
            if not cubes:
                continue
            srcs = srcs_for_field(target, prop_dir.name)
            try:
                row = src_df.loc[target]
                dist_kpc = float(row["dist_kpc"])
                target_ra = float(row["ra_deg"])
                target_dec = float(row["dec_deg"])
            except (KeyError, ValueError):
                dist_kpc, target_ra, target_dec = None, None, None
            for cube in cubes:
                n_files += 1
                orig, cut, status = maybe_cut(
                    cube, target, prop_dir.name, srcs, dist_kpc,
                    target_ra, target_dec, apply=args.apply
                )
                total_in += orig
                if cut is not None:
                    total_out += cut
                if status in ("wrote", "dryrun") or status.startswith("dryrun-"):
                    n_wrote += 1
                else:
                    n_skip += 1
                if n_files % 50 == 0:
                    print(f"  [{n_files}] {target}/{prop_dir.name}/"
                          f"{cube.name[:40]} {status} orig={orig/1e9:.2f}GB"
                          f"{f' cut={cut/1e9:.2f}GB' if cut else ''}",
                          flush=True)
    print(f"\n=== Summary ===")
    print(f"  Cubes scanned : {n_files}")
    print(f"  Cuts planned  : {n_wrote}")
    print(f"  Skipped       : {n_skip}")
    print(f"  Total in (GB) : {total_in/1e9:,.1f}")
    print(f"  Total out (GB): {total_out/1e9:,.1f}")
    if total_in:
        ratio = total_out / total_in
        print(f"  Reduction     : {(1-ratio)*100:.1f}%")


if __name__ == "__main__":
    main()
