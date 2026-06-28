"""Per-source detection-evidence audit.

For every analyzed (target, proposal, source) in analysis_products/, compile
a defensible detection / nondetection record for each species class plus
warning flags. Output:
  data/evidence_audit.csv

For each species (NaCl, KCl, H2O, RRL, SiO, SiS, SO):
  - lines_in_band: list of line names with in_band=True
  - n_5sig: count >= 5 sigma
  - top_line: line with max snr
  - top_snr, top_peak_K, top_peak_v_kms
  - stack_snr, stack_peak_K (if {NaCl,KCl,H2O}_stack.npz exists)
  - has_v0: whether any v=0 / ground-state line was measured
  - any_vexcited_5sig_with_v0_nondet: BIG red flag

Warning flags (column `flags`):
  H2O_NONDET_WITH_NACL_DET   - H2O lines in band but no 5sig + NaCl 5sig
  KCL_DET_NACL_NONDET        - KCl 5sig but NaCl not 5sig
  NACL_VEXC_ONLY             - NaCl v>=2 5sig but v=0 covered and <5sig
  KCL_VEXC_ONLY              - KCl v>=2 5sig but v=0 covered and <5sig
  SINGLE_LINE_DET            - exactly one in-band line drives the detection
  HIGH_CONFUSION             - aux f5sigma > 0.02 at this target+prop
"""
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
OUT_CSV = ROOT / "data/evidence_audit.csv"
AUX = ROOT / "data/data_summary_aux.csv"

SPECIES = {
    "NaCl": lambda L: L.startswith("NaCl_"),
    "KCl":  lambda L: L.startswith("KCl_") or L.startswith("KCl;"),
    "H2O":  lambda L: L.startswith("H2O"),
    "RRL":  lambda L: bool(re.match(r"^H\d+(alpha|beta|gamma|delta)$", L)),
    "SiO":  lambda L: L.startswith("SiO") or L.startswith("29SiO"),
    "SiS":  lambda L: L.startswith("SiS") or L.startswith("29SiS"),
    "SO":   lambda L: L.startswith("SO_") or L.startswith("SO2_"),
}

RE_VEXC = re.compile(r"_v(\d+)_")  # e.g. NaCl_v2_J17-16 -> '2'


def vibrational_state(line: str) -> int:
    m = RE_VEXC.search(line)
    if not m:
        return -1  # unknown
    return int(m.group(1))


def brightest_source_id(cont_csv):
    if not cont_csv.exists():
        return None
    try:
        df = pd.read_csv(cont_csv)
    except pd.errors.EmptyDataError:
        return None
    if df.empty or "peak_Jybeam" not in df.columns:
        return None
    return int(df.loc[df["peak_Jybeam"].idxmax(), "id"])


def stack_snr(prop_dir, sid, sp):
    if sp not in ("NaCl", "KCl"):
        return (np.nan, np.nan, 0)
    npz = prop_dir / f"source_{sid:02d}" / f"{sp}_stack.npz"
    if not npz.exists():
        return (np.nan, np.nan, 0)
    d = np.load(npz)
    if not {"spec", "sigma"} <= set(d.files):
        return (np.nan, np.nan, 0)
    spec = np.asarray(d["spec"])
    sigma = float(d["sigma"])
    nl = int(d["n_lines"]) if "n_lines" in d.files else 0
    if not np.isfinite(sigma) or sigma <= 0:
        return (np.nan, np.nan, nl)
    peak = float(np.nanmax(spec))
    return (peak / sigma, peak, nl)


def audit_one(prop_dir):
    cont = prop_dir / "continuum_sources.csv"
    bid = brightest_source_id(cont)
    if bid is None:
        return None
    meas_csv = prop_dir / "line_measurements.csv"
    if not meas_csv.exists():
        return None
    try:
        m = pd.read_csv(meas_csv)
    except pd.errors.EmptyDataError:
        return None
    if m.empty or "snr" not in m.columns:
        return None
    sub = m[m["source"] == bid].copy()
    target = prop_dir.parent.name
    proposal = prop_dir.name
    rec = {
        "target": target, "proposal": proposal, "src_id": bid,
    }
    flags = []
    for sp, pred in SPECIES.items():
        hits = sub[sub["line"].astype(str).apply(pred)].copy()
        # in-band: need not have measurement; treat presence in
        # line_measurements as 'in_band True'
        if hits.empty:
            for k in ("n_lines", "n_5sig", "top_line", "top_snr",
                      "top_peak_K", "top_peak_v",
                      "stack_snr", "stack_peak", "stack_nlines",
                      "vexc_only_det"):
                rec[f"{sp}_{k}"] = np.nan if k.startswith(("top_",
                                                              "stack_"))\
                    else 0 if k.startswith("n_") else ""
            rec[f"{sp}_decision"] = "na"
            continue
        hits_valid = hits[np.isfinite(hits["snr"])].copy()
        if hits_valid.empty:
            for k in ("n_lines", "n_5sig", "top_line", "top_snr",
                      "top_peak_K", "top_peak_v",
                      "stack_snr", "stack_peak", "stack_nlines",
                      "vexc_only_det"):
                rec[f"{sp}_{k}"] = np.nan if k.startswith(("top_", "stack_")) \
                    else (0 if k.startswith("n_") else "")
            rec[f"{sp}_decision"] = "na"
            continue
        n_5 = int((hits_valid["snr"] >= 5).sum())
        top = hits_valid.loc[hits_valid["snr"].idxmax()]
        hits = hits_valid
        s_snr, s_peak, s_n = stack_snr(prop_dir, bid, sp)

        # Vibrational excitation warning: 5sig line is v>=2 while v=0 covered
        # and not 5sig
        vexc_only = False
        if n_5 >= 1:
            states_5 = hits[hits["snr"] >= 5]["line"].apply(vibrational_state).tolist()
            states_all = hits["line"].apply(vibrational_state).tolist()
            high_v_5sig = any(v >= 2 for v in states_5)
            low_v_covered_nondet = (0 in states_all and
                                       not any(v == 0 for v in states_5))
            vexc_only = high_v_5sig and low_v_covered_nondet
            if vexc_only:
                flags.append(f"{sp}_VEXC_ONLY")

        # Single-line warning
        if n_5 == 1 and len(hits) >= 1:
            flags.append(f"{sp}_SINGLE_LINE_DET")

        rec[f"{sp}_n_lines"] = int(len(hits))
        rec[f"{sp}_n_5sig"] = n_5
        rec[f"{sp}_top_line"] = str(top["line"])
        rec[f"{sp}_top_snr"] = float(top["snr"])
        rec[f"{sp}_top_peak_K"] = float(top.get("peak_Kkms_or_unit", np.nan))
        rec[f"{sp}_top_peak_v"] = float(top.get("peak_v", np.nan))
        rec[f"{sp}_stack_snr"] = s_snr
        rec[f"{sp}_stack_peak"] = s_peak
        rec[f"{sp}_stack_nlines"] = s_n
        rec[f"{sp}_vexc_only_det"] = bool(vexc_only)
        # Decision: det if n_5sig >= 1 OR stack >= 4 sigma
        if n_5 >= 1 or (np.isfinite(s_snr) and s_snr >= 4.0):
            rec[f"{sp}_decision"] = "det"
        else:
            rec[f"{sp}_decision"] = "ul"

    # Cross-species warnings
    nacl_det = rec.get("NaCl_decision") == "det"
    kcl_det = rec.get("KCl_decision") == "det"
    h2o_det = rec.get("H2O_decision") == "det"
    h2o_covered = rec.get("H2O_decision", "na") != "na"
    nacl_covered = rec.get("NaCl_decision", "na") != "na"
    if nacl_det and h2o_covered and not h2o_det:
        flags.append("H2O_NONDET_WITH_NACL_DET")
    if kcl_det and nacl_covered and not nacl_det:
        flags.append("KCL_DET_NACL_NONDET")

    # Confusion flag from aux
    if AUX.exists():
        aux = pd.read_csv(AUX)
        a = aux[(aux["target"] == target) & (aux["proposal"] == proposal)]
        if not a.empty:
            f5 = float(a["f5sigma"].max())
            rec["f5sigma"] = f5
            if f5 > 0.02:
                flags.append("HIGH_CONFUSION")
        else:
            rec["f5sigma"] = np.nan
    else:
        rec["f5sigma"] = np.nan
    rec["flags"] = "|".join(flags)
    return rec


def main():
    rows = []
    for prop_dir in sorted(ANALYSIS.glob("*/2*")):
        r = audit_one(prop_dir)
        if r is None:
            continue
        rows.append(r)
    if not rows:
        print("nothing"); return
    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV} ({len(df)} rows)")
    # Pretty summary of flagged rows
    print("\n=== Sources with WARNING FLAGS ===")
    flagged = df[df["flags"] != ""]
    for _, r in flagged.iterrows():
        print(f"{r['target']:25s} {r['proposal']:18s} src{int(r['src_id']):02d}: {r['flags']}")


if __name__ == "__main__":
    main()
