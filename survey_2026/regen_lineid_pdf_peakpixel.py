"""Regenerate spectra/lineid_full.{png,pdf} for the new-detection targets,
using peak-pixel extraction at the brightest *on-target* mm source so the
spectrum actually contains the lines shown in the mom maps."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

sys.path.insert(0, str(Path(__file__).parent))
from analysis import diagnostics  # noqa: E402

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
TARGETS = ["G326.6618+00.5207", "G015.0357-00.6795",
            "G345.5043+00.3480", "IRAS17233-3606"]


def _catalog_coord(name):
    df = pd.read_csv(ROOT / "data/sources_L4_d2.csv")
    row = df[df["name"] == name]
    if row.empty:
        return None
    r = row.iloc[0]
    return float(r["ra_deg"]), float(r["dec_deg"])


def _on_target_peak(target, cat_coord):
    """Brightest mm source within 5\" of cat across this target's proposals.
    Returns (ra, dec) or None."""
    best = None
    best_peak = -1.0
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
            dra = (float(r["ra_deg"]) - cat_coord[0]) * np.cos(np.radians(cat_coord[1])) * 3600.0
            ddec = (float(r["dec_deg"]) - cat_coord[1]) * 3600.0
            sep = float(np.sqrt(dra * dra + ddec * ddec))
            if sep > 5.0:
                continue
            pk = float(r.get("peak_Jybeam", 0.0))
            if pk > best_peak:
                best_peak = pk
                best = (float(r["ra_deg"]), float(r["dec_deg"]))
    return best


def main():
    for tgt in TARGETS:
        rep = ANALYSIS / tgt / "report.json"
        if not rep.exists():
            print(f"  {tgt}: no report.json; skip")
            continue
        obj = json.loads(rep.read_text())
        cubes = obj.get("cube_paths", []) or []
        vcen = obj.get("vcen_kms", 0.0) or 0.0
        if not cubes:
            print(f"  {tgt}: no cube_paths; skip")
            continue
        cat = _catalog_coord(tgt) or _catalog_coord(tgt.replace("_", "-"))
        if cat is None:
            print(f"  {tgt}: no catalog coord; skip")
            continue
        peak = _on_target_peak(tgt, cat)
        if peak is None:
            print(f"  {tgt}: no on-target peak source found; falling back to catalog coord")
            ra, dec = cat
        else:
            ra, dec = peak
            print(f"  {tgt}: peak coord  ra={ra:.5f}  dec={dec:.5f}")
        coord = SkyCoord(ra * u.deg, dec * u.deg)
        outp = ANALYSIS / tgt / "spectra" / "lineid_full.png"
        try:
            diagnostics.make_lineid_full(
                cubes, str(outp), tgt, vcen,
                peak_coord=coord, mode="peak_pixel")
            print(f"  {tgt}: wrote lineid_full.png + .pdf (peak-pixel)")
        except (OSError, ValueError, RuntimeError) as e:
            print(f"  {tgt}: failed: {e}")


if __name__ == "__main__":
    main()
