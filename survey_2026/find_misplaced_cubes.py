"""Find ALMA cubes in uvdata/<proposal>/<target>/ that pointed at a DIFFERENT
target. For each off-target cube, look up the actual target name via SIMBAD
nearest-neighbor on the WCS pointing center, then either:

  --dry-run (default): report a move plan
  --move: actually shutil.move the file to uvdata/<proposal>/<simbad_name>/

A cube is considered "off-target" when its CRVAL1/CRVAL2 (pointing center)
is more than --tol-arcsec away from the brightest mm continuum source for
this target folder (or, when no analysis_products row exists, from the
target's L4_d2 coord).

For calibrator FITS files (quasar J<RA><DEC>, *_ph.*, *_bp.*, *_check.*),
the script moves them to uvdata/<proposal>/_calibrators/.
"""
import argparse
import re
import shutil
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
UVDIR = ROOT / "uvdata"
ANALYSIS = ROOT / "analysis_products"

CAL_PATTERNS = (
    re.compile(r"J\d{4}[-+]\d{4}"),
    re.compile(r"_ph\.", re.IGNORECASE),
    re.compile(r"_ampph\.", re.IGNORECASE),
    re.compile(r"_bp\.", re.IGNORECASE),
    re.compile(r"_check\.", re.IGNORECASE),
    re.compile(r"_pol\.", re.IGNORECASE),
)


def is_calibrator(filename: str) -> bool:
    return any(p.search(filename) for p in CAL_PATTERNS)


def brightest_source_coord(target):
    """Brightest mm continuum source across analyzed proposals for this
    target; fallback to L4_d2 coord."""
    tgt_dir = ANALYSIS / target
    best = None
    if tgt_dir.is_dir():
        for prop_dir in tgt_dir.glob("2*"):
            cont = prop_dir / "continuum_sources.csv"
            if not cont.exists():
                continue
            try:
                df = pd.read_csv(cont)
            except pd.errors.EmptyDataError:
                continue
            if df.empty:
                continue
            idx = int(df["peak_Jybeam"].idxmax())
            peak = float(df.loc[idx, "peak_Jybeam"])
            ra = float(df.loc[idx, "ra_deg"])
            dec = float(df.loc[idx, "dec_deg"])
            if best is None or peak > best[0]:
                best = (peak, ra, dec)
    if best is not None:
        return SkyCoord(best[1], best[2], unit="deg")
    # Fallback: L4_d2 coord
    src = pd.read_csv(ROOT / "data/sources_L4_d2.csv")
    r = src[src["name"] == target]
    if r.empty:
        return None
    return SkyCoord(float(r.iloc[0]["ra_deg"]), float(r.iloc[0]["dec_deg"]),
                      unit="deg")


def cube_pointing_coord(cube_path):
    h = fits.getheader(cube_path)
    return SkyCoord(float(h["CRVAL1"]), float(h["CRVAL2"]), unit="deg")


_simbad_cache: dict = {}


def simbad_lookup(coord, radius_arcsec=15):
    """Return the closest SIMBAD object name within radius_arcsec, or None."""
    key = (round(coord.ra.deg, 4), round(coord.dec.deg, 4), radius_arcsec)
    if key in _simbad_cache:
        return _simbad_cache[key]
    try:
        from astroquery.simbad import Simbad
        Simbad.TIMEOUT = 30
        tab = Simbad.query_region(coord, radius=radius_arcsec * u.arcsec)
    except (ConnectionError, TimeoutError, OSError, KeyError) as e:
        print(f"  SIMBAD fail: {e}")
        _simbad_cache[key] = None
        return None
    if tab is None or len(tab) == 0:
        _simbad_cache[key] = None
        return None
    # Pick brightest / closest. Simbad returns sorted by distance.
    name = str(tab[0]["main_id"]).strip()
    _simbad_cache[key] = name
    return name


def normalize_name(s):
    return re.sub(r"\s+", "", str(s)).replace("/", "_")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tol-arcsec", type=float, default=10.0,
                     help="Pointing-center tolerance for 'same target'.")
    ap.add_argument("--move", action="store_true",
                     help="Actually move files (default: dry-run).")
    ap.add_argument("--targets", nargs="*", default=None,
                     help="Limit to these L4_d2 target names.")
    args = ap.parse_args()

    src = pd.read_csv(ROOT / "data/sources_L4_d2.csv")
    target_coords = {r["name"]: SkyCoord(float(r["ra_deg"]),
                                            float(r["dec_deg"]),
                                            unit="deg")
                     for _, r in src.iterrows()}
    targets = args.targets or sorted(target_coords.keys())

    plan = []
    for target in targets:
        if target not in target_coords:
            continue
        ref = brightest_source_coord(target) or target_coords[target]
        for prop_dir in sorted(UVDIR.iterdir()):
            sub = prop_dir / target
            if not sub.is_dir():
                continue
            for f in sub.iterdir():
                if not f.name.endswith(".pbcor.fits"):
                    continue
                if is_calibrator(f.name):
                    new_dir = prop_dir / "_calibrators"
                    plan.append((f, new_dir / f.name, "calibrator"))
                    continue
                try:
                    pc = cube_pointing_coord(f)
                except (KeyError, OSError):
                    continue
                sep = ref.separation(pc).arcsec
                if sep < args.tol_arcsec:
                    continue
                # Check if any other L4_d2 target is closer
                seps_all = {n: c.separation(pc).arcsec
                            for n, c in target_coords.items()}
                nearest = min(seps_all.items(), key=lambda kv: kv[1])
                if nearest[1] < args.tol_arcsec:
                    new_dir = prop_dir / nearest[0]
                    plan.append((f, new_dir / f.name,
                                  f"matches L4_d2 {nearest[0]} ({nearest[1]:.1f}\")"))
                    continue
                # SIMBAD lookup as last resort
                name = simbad_lookup(pc)
                if name:
                    safe = normalize_name(name)
                    new_dir = prop_dir / safe
                    plan.append((f, new_dir / f.name,
                                  f"SIMBAD: {name} (sep {sep:.0f}\")"))
                else:
                    plan.append((f, None, f"unknown target, off by {sep:.0f}\""))

    print(f"\n=== Move plan ({len(plan)} files) ===")
    for src_p, dest_p, note in plan:
        print(f"  {src_p.relative_to(UVDIR)}")
        print(f"    -> {dest_p.relative_to(UVDIR) if dest_p else 'KEEP'}  [{note}]")
    if args.move:
        n_moved = 0
        for src_p, dest_p, _ in plan:
            if dest_p is None:
                continue
            dest_p.parent.mkdir(parents=True, exist_ok=True)
            if dest_p.exists():
                print(f"  exists, skip: {dest_p}")
                continue
            shutil.move(str(src_p), str(dest_p))
            n_moved += 1
        print(f"\nMoved {n_moved} files.")
    else:
        print("\n(dry-run: pass --move to actually move)")


if __name__ == "__main__":
    main()
