"""Apply strict-vet criteria to every target in analysis_products and
contrast against build_data_summary's T4 cell flags.

For each on-target source (within 5" of catalog), examine every salt/H2O
spec.npz and compute the same flags vet_new_detections uses (peak SNR,
trough_clean, v_ok). Per-species verdict (REAL / SUSPECT / NON-DETECTION)
is then compared against what the current T4 cell says (y/n) so we can
see which targets would be promoted or demoted.

Writes data/all_detection_vet.csv and prints a per-target diff summary.
"""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
OUT_CSV = ROOT / "data/all_detection_vet.csv"

SPECIES_LINES = {
    "NaCl": "NaCl_*.spec.npz",
    "KCl":  "KCl_*.spec.npz",
    "H2O":  "H2O_*.spec.npz",
}


def _vlsr_for(target):
    for f in ("data/vlsr_from_data.json", "data/vlsr_from_literature.json"):
        p = ROOT / f
        if not p.exists():
            continue
        try:
            obj = json.loads(p.read_text())
        except json.JSONDecodeError:
            continue
        rec = (obj.get(target) or obj.get(target.replace("-", "_"))
               or obj.get(target.replace("_", "-")))
        if rec and rec.get("v_LSR_kms") is not None:
            return float(rec["v_LSR_kms"])
    return None


def _line_label(name):
    parts = name.split("_")
    vn = None
    for p in parts[1:]:
        if p.startswith("v") and len(p) > 1 and p[1:].isdigit():
            vn = int(p[1:])
    if name.startswith("H2O_") and vn is None:
        vn = 0
    return vn


def _on_target_srcs(target, ra_cat, dec_cat):
    """Yield (proposal_dir, src_id) for sources within 5" of catalog."""
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
                yield pdir, int(r["id"])


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


def main():
    src_df = pd.read_csv(ROOT / "data/sources_L4_d2.csv")
    cat = {r["name"]: (float(r["ra_deg"]), float(r["dec_deg"]))
           for _, r in src_df.iterrows()}
    rows = []
    for tdir in sorted(ANALYSIS.iterdir()):
        if not tdir.is_dir() or tdir.name.startswith("_"):
            continue
        target = tdir.name
        coord = (cat.get(target) or cat.get(target.replace("_", "-"))
                  or cat.get(target.replace("-", "_")))
        if coord is None:
            continue
        ra_cat, dec_cat = coord
        vlsr = _vlsr_for(target)
        for pdir, sid in _on_target_srcs(target, ra_cat, dec_cat):
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
                    vn = _line_label(npz.stem)
                    dv = (st["peak_v"] - vlsr) if vlsr is not None else 0.0
                    v_ok = abs(dv) <= 10.0
                    trough_clean = (tsn < 0.5 * psn) and (tsn < 5.0)
                    line_name = npz.stem[:-5] if npz.stem.endswith(".spec") else npz.stem
                    rows.append(dict(
                        target=target, proposal=pdir.name, src_id=sid,
                        species=species, line=line_name, v=vn,
                        peak_K=st["peak"], trough_K=st["trough"],
                        sigma_K=st["sigma"], snr=psn, trough_snr=tsn,
                        peak_v=st["peak_v"], dv_from_vlsr=dv,
                        v_ok=v_ok, trough_clean=trough_clean, is_v0=(vn == 0),
                        passes_strict=(psn >= 5.0 and trough_clean and v_ok),
                    ))
    if not rows:
        print("no rows"); return
    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV}  ({len(df)} rows)")
    # Per-target diff vs T4 inferred from passes
    # Verdict = REAL if v=0 pass OR >=2 passes; SUSPECT if 1 pass; ND otherwise.
    print()
    print(f"{'target':<25} {'species':<5} {'verdict':<14} {'n_pass/n_lines':<14} {'v0?':<5}")
    for (tgt, sp), g in df.groupby(["target", "species"]):
        # collapse over all proposals/sources for that target+species
        np_ = int(g["passes_strict"].sum())
        v0_det = bool(((g["is_v0"]) & (g["passes_strict"])).any())
        verdict = ("REAL" if (v0_det or np_ >= 2)
                   else "SUSPECT" if np_ >= 1 else "ND")
        if verdict != "ND":  # only print interesting rows
            print(f"{tgt:<25} {sp:<5} {verdict:<14} {np_}/{len(g):<13} {v0_det}")


if __name__ == "__main__":
    main()
