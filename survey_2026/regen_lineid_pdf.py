"""Regenerate spectra/lineid_full.{png,pdf} for a fixed list of targets."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from analysis import diagnostics  # noqa: E402

ANALYSIS = Path("/orange/adamginsburg/salt/survey_2026/analysis_products")

TARGETS = [
    "G326.6618+00.5207",
    "G015.0357-00.6795",
    "G345.5043+00.3480",
    "IRAS17233-3606",
]


def main():
    for tgt in TARGETS:
        src_dir = ANALYSIS / tgt
        rep_path = src_dir / "report.json"
        if not rep_path.exists():
            print(f"  {tgt}: no report.json; skip")
            continue
        rep = json.loads(rep_path.read_text())
        cube_paths = rep.get("cube_paths", []) or []
        vcen = rep.get("vcen_kms", 0.0) or 0.0
        if not cube_paths:
            print(f"  {tgt}: no cube_paths in report.json; skip")
            continue
        outp = src_dir / "spectra" / "lineid_full.png"
        try:
            diagnostics.make_lineid_full(cube_paths, str(outp), tgt, vcen)
            print(f"  {tgt}: wrote {outp.name} + .pdf")
        except (OSError, ValueError, RuntimeError) as e:
            print(f"  {tgt}: failed: {e}")


if __name__ == "__main__":
    main()
