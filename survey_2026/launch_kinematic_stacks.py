"""Launch kinematic stacks for every (target, src) with at least one
probable, high-SNR (snr>=8) detection from the per_source_linelist.csv.

For each (target), pick the highest-SNR probable line as the velocity
guide and spawn `python build_kinematic_stack.py --target ...`. One
sub-process per target; runs sequentially.
"""
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
LINELIST = ROOT / "data/per_source_linelist.csv"
VLSR_LIT = ROOT / "data/vlsr_from_literature.json"
VLSR_DATA = ROOT / "data/vlsr_from_data.json"


def vlsr_for(target):
    for p in (VLSR_DATA, VLSR_LIT):
        if not p.exists():
            continue
        try:
            d = json.loads(p.read_text())
            if target in d and d[target].get("v_LSR_kms") is not None:
                return float(d[target]["v_LSR_kms"])
        except (json.JSONDecodeError, AttributeError):
            pass
    return None


def main():
    df = pd.read_csv(LINELIST)
    df = df[(df["probable"]) & (df["snr"] >= 8.0)]
    if df.empty:
        print("no candidates")
        return
    # Best guide per (target, proposal, src) — top SNR row
    df = df.sort_values("snr", ascending=False)
    seen = set()
    plans = []
    for _, r in df.iterrows():
        key = (r["target"], r["proposal"], int(r["src"]))
        if key in seen:
            continue
        seen.add(key)
        plans.append(r)
    print(f"{len(plans)} kinematic stacks to launch")
    for r in plans:
        target = str(r["target"])
        prop = str(r["proposal"])
        src = int(r["src"])
        guide = str(r["line"])
        rest = float(r["rest_GHz"])
        vlsr = vlsr_for(target)
        if vlsr is None:
            # Fall back to the guide line's measured peak velocity
            try:
                vlsr = float(r["peak_v"])
                print(f"  {target}: no lit vLSR, using line peak_v="
                      f"{vlsr:+.1f} km/s")
            except (KeyError, TypeError, ValueError):
                print(f"  SKIP {target}: no vLSR")
                continue
        cmd = [sys.executable, "build_kinematic_stack.py",
               "--target", target, "--proposal", prop,
               "--guide-line", guide, "--guide-rest-GHz", f"{rest:.6f}",
               "--vlsr", f"{vlsr:.2f}", "--source-id", str(src)]
        print(f"\n=== {target}/{prop} src{src:02d} guide={guide} "
              f"vlsr={vlsr:+.1f} ===", flush=True)
        try:
            subprocess.run(cmd, check=False, timeout=600,
                              cwd=str(ROOT))
        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT after 600 s")


if __name__ == "__main__":
    main()
