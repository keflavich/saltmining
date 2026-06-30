"""Vet CARTA snippets against source_conflict_audit.csv.

For each saltsearch2026 snippet, parse the openFile/appendFile/importRegion
lines, drop any line referencing an off-target (target, proposal) pair as
flagged in the audit. Writes back a cleaned snippet. Also reports paths
that point to files that no longer exist on disk.
"""
import json
import re
import warnings
from pathlib import Path

import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
SNIPPET_DIR = Path.home() / ".carta/config/snippets"
AUDIT_CSV = ROOT / "data/source_conflict_audit.csv"
CARTA_ROOT = Path("/orange/adamginsburg")  # snippets store relative paths

PATH_RE = re.compile(r'"([^"]+)"')


def main():
    if not AUDIT_CSV.exists():
        print(f"audit csv missing: {AUDIT_CSV}; run audit_source_conflicts.py first")
        return
    audit = pd.read_csv(AUDIT_CSV)
    # Set of (target, proposal) pairs that are entirely off-target.
    off_target = set()
    for (t, p), grp in audit.groupby(["target", "proposal"]):
        if (grp["all_off_target"]).all() and not grp["on_target"].any():
            off_target.add((t, p))
    print(f"off-target proposals: {len(off_target)}")

    fixed = 0
    for snip in sorted(SNIPPET_DIR.glob("*_2026.json")):
        try:
            obj = json.loads(snip.read_text())
        except json.JSONDecodeError:
            continue
        if "saltsearch2026" not in (obj.get("categories") or []):
            continue
        code = obj.get("code", "") or ""
        lines = code.split("\n")
        kept = []
        dropped = []
        for line in lines:
            paths = PATH_RE.findall(line)
            keep = True
            for path in paths:
                # Translate carta-relative path back to absolute
                abs_path = str(CARTA_ROOT / path.lstrip("/"))
                # Check off-target: path must contain /<proposal>/<target>/
                for (t, p) in off_target:
                    # uvdata layout: /<prop>/<target>/
                    # analysis_products layout: /<target>/<prop>/
                    if (f"/{p}/{t}/" in abs_path
                        or f"/{t}/{p}/" in abs_path):
                        keep = False
                        dropped.append((line, t, p))
                        break
                # Also drop if file no longer exists on disk (e.g.,
                # quarantined into _misidentified/).
                ap = Path(abs_path)
                if ap.suffix in (".fits", ".reg") and not ap.exists():
                    keep = False
                    dropped.append((line, "missing", ""))
                if not keep:
                    break
            if keep:
                kept.append(line)
        if dropped:
            print(f"\n{snip.name}: dropping {len(dropped)} lines")
            for line, t, p in dropped[:3]:
                print(f"  - {t}/{p}: {line[:80]}")
            new_code = "\n".join(kept)
            if new_code.strip():
                obj["code"] = new_code
                snip.write_text(json.dumps(obj, indent=4))
                fixed += 1
            else:
                # All lines dropped -> remove the snippet outright
                snip.unlink()
                print(f"  removed empty snippet")
                fixed += 1
    print(f"\nupdated {fixed} snippet(s)")


if __name__ == "__main__":
    main()
