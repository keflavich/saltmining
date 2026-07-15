"""Run build_pv_diagrams over every (target, proposal) with a strict-vet
REAL detection in data/all_detection_vet.csv. Skips targets without a
guide-line mom1 already on disk."""
import json
import subprocess
import warnings
from pathlib import Path

import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"


def _vlsr(target):
    for f in ("data/vlsr_from_data.json", "data/vlsr_from_literature.json"):
        p = ROOT / f
        if not p.exists():
            continue
        try:
            d = json.loads(p.read_text())
        except json.JSONDecodeError:
            continue
        rec = (d.get(target) or d.get(target.replace("-", "_"))
               or d.get(target.replace("_", "-")))
        if rec and rec.get("v_LSR_kms") is not None:
            return float(rec["v_LSR_kms"])
    # fall back to the survey-wide obs_params vsrc (covers all targets)
    import pandas as pd
    obs = pd.read_csv(ROOT / "data/obs_params.csv")
    row = obs[obs["name"] == target]
    if len(row) and pd.notna(row.iloc[0]["vsrc_kms"]):
        return float(row.iloc[0]["vsrc_kms"])
    print(f"WARNING: no vlsr for {target}; using 0.0")
    return 0.0


def main():
    vet = pd.read_csv(ROOT / "data/all_detection_vet.csv")
    real = vet[vet["passes_strict"]]
    # Aggregate (target, proposal) pairs that have at least one REAL line.
    pairs = real.groupby(["target", "proposal"]).size().reset_index(name="n")
    pairs = pairs.sort_values("n", ascending=False)
    print(f"{len(pairs)} (target, proposal) pairs with REAL passes")
    for _, r in pairs.iterrows():
        target = str(r["target"]); proposal = str(r["proposal"])
        vlsr = _vlsr(target)
        # Look for a usable guide line: prefer H2O 5_15, then NaCl_v0, then any
        guide = None
        sdir_glob = list((ANALYSIS / target / proposal).glob("source_*"))
        for s in sdir_glob:
            for pref in ("H2O_5_15-4_22_232", "NaCl_v0_J18-17",
                         "H2O_232p9769", "H30alpha"):
                if (s / f"{s.name}_{pref}_mom1.fits").exists():
                    guide = pref
                    break
            if guide:
                break
        if guide is None:
            print(f"  {target} {proposal}: no usable guide-line mom1; skip")
            continue
        print(f"\n=== {target}  {proposal}  vlsr={vlsr:+.1f}  guide={guide}")
        cmd = ["python", "-u", "build_pv_diagrams.py",
               "--target", target,
               "--proposal", proposal,
               "--vlsr", str(vlsr),
               "--guide-line", guide]
        try:
            subprocess.run(cmd, cwd=str(ROOT), timeout=900, check=False)
        except subprocess.TimeoutExpired:
            print(f"  ! timeout")


if __name__ == "__main__":
    main()
