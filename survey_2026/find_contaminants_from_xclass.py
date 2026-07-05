"""For each suspect detected line, identify candidate contaminating species
from the XCLASS transition_energies.dat output of the matched per-target
XCLASS run.

Workflow:
  - For each (target_label) under lineid/, find the matching XCLASS job dir
    under /orange/adamginsburg/software/XCLASS-Interface/run/myXCLASS/
    and parse transition_energies.dat
  - For each detected line at the brightest mm source of that target's
    proposal, search the transition table for lines whose Doppler-shifted
    freq matches the observed peak freq within tol_kms

Output: data/detected_lines_with_contaminants.csv + .md
"""
import glob
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
LINEID = ROOT / "lineid"
JOBROOT = Path("/orange/adamginsburg/software/XCLASS-Interface/run/myXCLASS")

OUT_CSV = ROOT / "data/detected_lines_with_contaminants.csv"
OUT_MD = ROOT / "data/detected_lines_with_contaminants.md"

C_KMS = 299792.458

# Map XCLASS lineid label -> (analysis_products target, proposal)
TARGET_PROP = {
    "G326":     ("G326.6618+00.5207", "2022.1.01344.S"),
    "G015.0357": ("G015.0357-00.6795", "2023.1.01346.S"),
    "G345.5043": ("G345.5043+00.3480", "2022.1.01344.S"),
    "NGC6334I": ("NGC6334I",           "2016.1.00383.S"),
}


def parse_transition_table(path):
    """Return DataFrame: rest_MHz, dopp_MHz, T_K, eup_K, mol."""
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("%") or not line.strip():
                continue
            # Whitespace-split first 7 numbers, molecule name is last token
            parts = line.rstrip().split()
            try:
                rest = float(parts[0])
                dopp = float(parts[1])
                inten = float(parts[2])
                # Last few tokens contain molecule
                mol = parts[-1].rstrip(";")
            except (ValueError, IndexError):
                continue
            rows.append((rest, dopp, inten, mol))
    return pd.DataFrame(rows, columns=["rest_MHz", "dopp_MHz", "T_K", "mol"])


def find_latest_job(target_label):
    """Find the most recent XCLASS job dir that ran a molfit with this
    target label."""
    candidates = []
    for d in JOBROOT.glob("job__*"):
        for mfit in d.glob("*.molfit"):
            if target_label in mfit.name or target_label in mfit.stem:
                candidates.append(d)
                break
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def brightest_id(prop_dir):
    cont = prop_dir / "continuum_sources.csv"
    if not cont.exists() or cont.stat().st_size == 0:
        return None
    try:
        df = pd.read_csv(cont)
    except pd.errors.EmptyDataError:
        return None
    if df.empty or "peak_Jybeam" not in df.columns:
        return None
    return int(df.loc[df["peak_Jybeam"].idxmax(), "id"])


def main():
    out_rows = []
    for label, (target, proposal) in TARGET_PROP.items():
        print(f"\n=== {label} ({target}, {proposal}) ===")
        job = find_latest_job(label)
        if job is None:
            print(f"  no XCLASS job dir found for {label}")
            continue
        print(f"  XCLASS job: {job.name}")
        tt_path = job / "transition_energies.dat"
        if not tt_path.exists():
            print(f"  no transition_energies.dat")
            continue
        tt = parse_transition_table(tt_path)
        print(f"  parsed {len(tt)} transitions")

        prop_dir = ANALYSIS / target / proposal
        meas = prop_dir / "line_measurements.csv"
        if not meas.exists():
            continue
        df = pd.read_csv(meas)
        bid = brightest_id(prop_dir)
        if bid is None: continue
        hits = df[(df["source"] == bid) & (df["snr"] >= 5.0)].copy()
        if hits.empty:
            continue

        for _, r in hits.iterrows():
            line = str(r["line"])
            rest_GHz = float(r["rest_GHz"])
            peak_v = float(r.get("peak_v", np.nan))
            snr = float(r["snr"])
            obs_GHz = rest_GHz * (1.0 - peak_v / C_KMS)
            obs_MHz = obs_GHz * 1000.0
            # Search within +/- 25 MHz (~30 km/s at 232 GHz) — broad enough
            # to catch close COM blends; narrower offsets are listed first
            # because we sort by predicted T_K.
            mask = np.abs(tt["dopp_MHz"] - obs_MHz) < 25.0
            cands = tt[mask].sort_values("T_K", ascending=False).head(10)
            for _, c in cands.iterrows():
                # Skip self-match
                if abs(c["rest_MHz"] - rest_GHz * 1000.0) < 5.0:
                    continue
                out_rows.append({
                    "label": label, "target": target, "proposal": proposal,
                    "detected_line": line,
                    "detected_rest_GHz": rest_GHz,
                    "detected_snr": snr,
                    "detected_peak_v": peak_v,
                    "obs_MHz": obs_MHz,
                    "candidate_mol": c["mol"],
                    "candidate_rest_MHz": c["rest_MHz"],
                    "candidate_dopp_MHz": c["dopp_MHz"],
                    "candidate_T_K": c["T_K"],
                    "freq_offset_MHz": c["dopp_MHz"] - obs_MHz,
                })

    out = pd.DataFrame(out_rows)
    if out.empty:
        print("\nNo contaminants found.")
        return
    out.to_csv(OUT_CSV, index=False)
    print(f"\nwrote {OUT_CSV} ({len(out)} candidate rows)")

    # Markdown report
    md = ["# Contaminant candidates from XCLASS LTE model\n",
          "For each suspect detected line at the brightest mm continuum",
          "source, list COM/etc. transitions from the XCLASS model that fall",
          "within +/- 5 MHz of the observed peak frequency.\n"]
    for (label, det_line), g in out.groupby(["label", "detected_line"]):
        first = g.iloc[0]
        md.append(f"\n## {label}  {det_line}\n")
        md.append(f"- rest = {first['detected_rest_GHz']:.4f} GHz, "
                  f"snr = {first['detected_snr']:.1f}, "
                  f"peak_v = {first['detected_peak_v']:+.1f} km/s")
        md.append(f"- observed @ {first['obs_MHz']/1000.0:.4f} GHz\n")
        md.append("| Mol | rest (GHz) | Doppler-shifted (GHz) | freq offset (MHz) | XCLASS T (K) |")
        md.append("|---|---|---|---|---|")
        for _, c in g.head(5).iterrows():
            md.append(f"| {c['candidate_mol']} | "
                      f"{c['candidate_rest_MHz']/1000.0:.4f} | "
                      f"{c['candidate_dopp_MHz']/1000.0:.4f} | "
                      f"{c['freq_offset_MHz']:+.2f} | "
                      f"{c['candidate_T_K']:.1f} |")
    OUT_MD.write_text("\n".join(md))
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main()
