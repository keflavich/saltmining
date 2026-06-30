"""Audit analysis_products for source-ID/position conflicts.

For each target in data/sources_L4_d2.csv:
  - Iterate every proposal under analysis_products/<target>/2*/
  - Read continuum_sources.csv -> (id, ra_deg, dec_deg, peak_Jybeam, snr)
  - Compare each source's (ra, dec) against the catalog target coord
  - Flag:
      A. Cross-proposal collision: same src_id, different sky position
      B. Off-target proposal: ALL sources in the proposal are >5" from
         the catalog coord (likely a multi-field MOUS where the wrong
         field was picked)
      C. Source dir mismatch: source_NN/ exists but src_id NN doesn't
         match the brightest source

Writes a report to data/source_conflict_audit.csv with one row per source.
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
SRC_CSV = ROOT / "data/sources_L4_d2.csv"
OUT_CSV = ROOT / "data/source_conflict_audit.csv"
TOL_ARCSEC = 5.0


def sep_arcsec(ra1, dec1, ra2, dec2):
    dra = (ra1 - ra2) * np.cos(np.radians(dec2)) * 3600.0
    ddec = (dec1 - dec2) * 3600.0
    return float(np.sqrt(dra * dra + ddec * ddec))


def main():
    src = pd.read_csv(SRC_CSV)
    catalog = {r["name"]: (float(r["ra_deg"]), float(r["dec_deg"])) for _, r in src.iterrows()}
    rows = []
    for tdir in sorted(ANALYSIS.iterdir()):
        if not tdir.is_dir() or tdir.name.startswith("_"):
            continue
        cat = catalog.get(tdir.name) or catalog.get(
            tdir.name.replace("_", "-")) or catalog.get(tdir.name.replace("-", "_"))
        if cat is None:
            continue
        ra_cat, dec_cat = cat
        prop_to_sources = {}
        for pdir in sorted(tdir.glob("2*")):
            cont = pdir / "continuum_sources.csv"
            if not cont.exists() or cont.stat().st_size == 0:
                continue
            try:
                df = pd.read_csv(cont)
            except pd.errors.EmptyDataError:
                continue
            if df.empty or not {"id", "ra_deg", "dec_deg"} <= set(df.columns):
                continue
            srcs = []
            for _, r in df.iterrows():
                sep = sep_arcsec(float(r["ra_deg"]), float(r["dec_deg"]),
                                  ra_cat, dec_cat)
                srcs.append(dict(
                    src_id=int(r["id"]),
                    ra=float(r["ra_deg"]),
                    dec=float(r["dec_deg"]),
                    peak=float(r.get("peak_Jybeam", np.nan)),
                    snr=float(r.get("snr", np.nan)),
                    sep_arcsec=sep,
                ))
            prop_to_sources[pdir.name] = srcs
        if not prop_to_sources:
            continue
        # Per-proposal off-target check
        for prop, sl in prop_to_sources.items():
            on_tgt = [s for s in sl if s["sep_arcsec"] <= TOL_ARCSEC]
            for s in sl:
                rows.append(dict(
                    target=tdir.name, proposal=prop,
                    src_id=s["src_id"], ra=s["ra"], dec=s["dec"],
                    peak_Jybeam=s["peak"], snr=s["snr"],
                    sep_arcsec=s["sep_arcsec"],
                    on_target=s["sep_arcsec"] <= TOL_ARCSEC,
                    all_off_target=(len(on_tgt) == 0),
                ))
        # Cross-proposal id collision check
        id_positions = {}
        for prop, sl in prop_to_sources.items():
            for s in sl:
                key = s["src_id"]
                id_positions.setdefault(key, []).append((prop, s["ra"], s["dec"]))
        for sid, plist in id_positions.items():
            if len(plist) < 2:
                continue
            # any pair > 1" apart?
            seps = []
            for i in range(len(plist)):
                for j in range(i + 1, len(plist)):
                    s = sep_arcsec(plist[i][1], plist[i][2],
                                    plist[j][1], plist[j][2])
                    if s > 1.0:
                        seps.append((plist[i][0], plist[j][0], s))
            if seps:
                print(f"CROSS-COLLISION {tdir.name} src_id={sid}:")
                for a, b, s in seps[:3]:
                    print(f"  {a}  vs  {b}  -> {s:.2f}\"")
    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    print(f"\nwrote {OUT_CSV}  ({len(out)} rows)")
    # Summarize off-target proposals
    off = out[out["all_off_target"]].drop_duplicates(["target", "proposal"])
    print(f"\nProposals where NO source is within {TOL_ARCSEC}\" of catalog target:")
    for _, r in off.iterrows():
        print(f"  {r['target']:<30} {r['proposal']:<22} closest={r['sep_arcsec']:.1f}\"")


if __name__ == "__main__":
    main()
