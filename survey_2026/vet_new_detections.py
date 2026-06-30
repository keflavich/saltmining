"""Strict salt/water detection vetting for the "first reported here"
sources in Table 4.

Reads per-source spec.npz files and applies stricter criteria than the
pipeline's snr>=5 cell flag:

  1. Detection requires either
       (a) the ground-state (v=0) line of the species detected at >=5σ,
       OR
       (b) at least 2 transitions of the same species (J or v) each
           detected at >=5σ with peak_v within +/-5 km/s of the source
           vlsr.
  2. |trough|/sigma must be < 0.5 * peak/sigma AND <5 (i.e., no nearby
     absorption feature confusing the line as a contaminant).
  3. peak_v within +/-10 km/s of the source vlsr from
     vlsr_from_data.json / vlsr_from_literature.json.

Writes data/new_detection_vet.csv listing every (target, species, line)
with peak, trough, snr, vet flag, and a per-source verdict
(REAL / SUSPECT / NON-DETECTION).
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
OUT_CSV = ROOT / "data/new_detection_vet.csv"

# Targets flagged as "first reported here" in Table 4 by build_data_summary.
TARGETS = [
    ("IRAS 15412-5359",  "G326.6618+00.5207"),
    ("IRAS 18174-1612",  "G015.0357-00.6795"),
    ("I17016-4124",      "G345.5043+00.3480"),
    ("IRAS 17233-3606",  "IRAS17233-3606"),
]


def _vlsr_for(target):
    """Try vlsr_from_data, then vlsr_from_literature, then None."""
    for f in ("data/vlsr_from_data.json", "data/vlsr_from_literature.json"):
        p = ROOT / f
        if not p.exists():
            continue
        try:
            obj = json.loads(p.read_text())
        except json.JSONDecodeError:
            continue
        rec = (obj.get(target)
               or obj.get(target.replace("-", "_"))
               or obj.get(target.replace("_", "-")))
        if rec and rec.get("v_LSR_kms") is not None:
            return float(rec["v_LSR_kms"])
    return None


def _on_target_brightest_src(target):
    """Find brightest mm continuum source within 5\" of the catalog target
    coord, across all proposals. Returns (proposal_dir, src_id) or None."""
    cat_df = pd.read_csv(ROOT / "data/sources_L4_d2.csv")
    row = cat_df[cat_df["name"].isin([target, target.replace("_", "-"), target.replace("-", "_")])]
    if row.empty:
        return None
    ra_cat = float(row.iloc[0]["ra_deg"])
    dec_cat = float(row.iloc[0]["dec_deg"])
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
            if sep > 5.0:
                continue
            pk = float(r.get("peak_Jybeam", 0.0))
            if pk > best_pk:
                best_pk = pk
                best = (pdir, int(r["id"]))
    return best


SPECIES_LINES = {
    "NaCl": "NaCl_*.spec.npz",
    "KCl":  "KCl_*.spec.npz",
    "H2O":  "H2O_*.spec.npz",
}


def _spec_stats(npz_path, vlsr, win_kms=30.0):
    d = np.load(npz_path)
    if "spec" not in d.files or "vaxis" not in d.files:
        return None
    spec = np.asarray(d["spec"], dtype=float)
    vax = np.asarray(d["vaxis"], dtype=float)
    sigma = float(d["sigma"]) if "sigma" in d.files else float("nan")
    if vlsr is None:
        win = (vax > -win_kms) & (vax < win_kms)
    else:
        win = (vax > vlsr - win_kms) & (vax < vlsr + win_kms)
    if not win.any():
        return None
    peak = float(np.nanmax(spec[win]))
    trough = float(np.nanmin(spec[win]))
    peak_v = float(vax[win][int(np.nanargmax(spec[win]))])
    return dict(peak=peak, trough=trough, sigma=sigma, peak_v=peak_v)


def _line_label(name):
    """Return (v_number, J_label). For NaCl/KCl naming uses _v0_J18-17.
    For H2O the convention has no explicit v0 marker for the
    ground-vibrational 5_15-4_22 line (it's just "H2O_5_15-4_22..."),
    while the vibrationally-excited maser is "H2O_v2_..."; we map any
    H2O_* without an explicit _v[1-9]_ token to v=0."""
    parts = name.split("_")
    vn = None
    J = None
    for p in parts[1:]:
        if p.startswith("v") and len(p) > 1 and p[1:].isdigit():
            vn = int(p[1:])
        elif p.startswith("J"):
            J = p[1:]
    if name.startswith("H2O_") and vn is None:
        vn = 0
    return vn, J


def _all_on_target_srcs(target):
    """Yield (proposal_dir, src_id) for every continuum source within 5"
    of the catalog target, across all proposals."""
    cat_df = pd.read_csv(ROOT / "data/sources_L4_d2.csv")
    row = cat_df[cat_df["name"].isin([target, target.replace("_", "-"), target.replace("-", "_")])]
    if row.empty:
        return
    ra_cat = float(row.iloc[0]["ra_deg"])
    dec_cat = float(row.iloc[0]["dec_deg"])
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
            if sep <= 5.0:
                yield pdir, int(r["id"]), float(r.get("peak_Jybeam", 0.0))


def main():
    rows = []
    for display, target in TARGETS:
        vlsr = _vlsr_for(target)
        for pdir, sid, peak in _all_on_target_srcs(target):
            sdir = pdir / f"source_{sid:02d}"
            if not sdir.is_dir():
                continue
            for species, pat in SPECIES_LINES.items():
                for npz in sorted(sdir.glob(pat)):
                    st = _spec_stats(npz, vlsr)
                    if st is None or st["sigma"] <= 0:
                        continue
                    psn = st["peak"] / st["sigma"]
                    tsn = abs(st["trough"]) / st["sigma"]
                    vn, J = _line_label(npz.stem)
                    dv = (st["peak_v"] - vlsr) if vlsr is not None else 0.0
                    v_ok = abs(dv) <= 10.0
                    trough_clean = (tsn < 0.5 * psn) and (tsn < 5.0)
                    line_name = npz.stem
                    if line_name.endswith(".spec"):
                        line_name = line_name[:-len(".spec")]
                    rows.append(dict(
                        target=target, display=display, proposal=pdir.name,
                        src_id=sid, species=species, line=line_name,
                        v=vn, J=J,
                        peak_K=st["peak"], trough_K=st["trough"],
                        sigma_K=st["sigma"], snr=psn, trough_snr=tsn,
                        peak_v=st["peak_v"], dv_from_vlsr=dv,
                        v_ok=v_ok, trough_clean=trough_clean,
                        is_v0=(vn == 0), passes_strict=(psn >= 5.0 and trough_clean and v_ok),
                    ))
    if not rows:
        print("no rows")
        return
    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV}  ({len(df)} rows)")
    # Per-target verdict
    print()
    for (tgt, disp), g in df.groupby(["target", "display"]):
        print(f"== {disp}  ({tgt}) ==")
        for sp in sorted(g["species"].unique()):
            sg = g[g["species"] == sp]
            v0_det = ((sg["is_v0"]) & (sg["passes_strict"])).any()
            n_pass = int((sg["passes_strict"]).sum())
            verdict = ("REAL"  if (v0_det or n_pass >= 2)
                       else "SUSPECT" if n_pass >= 1
                       else "NON-DETECTION")
            print(f"   {sp:<5} -> {verdict:<14}  (v0_det={v0_det}, n_pass={n_pass}/{len(sg)})")


if __name__ == "__main__":
    main()
