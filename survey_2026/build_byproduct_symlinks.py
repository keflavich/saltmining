"""Build a parallel `by_type/` symlink tree under analysis_products/ that
groups generated plots by product class instead of by (target, proposal).

Outputs:

  analysis_products/by_type/overview_continuum/<target>_<proposal>.png
  analysis_products/by_type/spectrum_full/<target>_<proposal>_src<NN>[_p<N>].png
  analysis_products/by_type/kinematic_stack/<target>_<proposal>_<guide>.png
  analysis_products/by_type/salt_stack/<target>_<proposal>_src<NN>.png
  analysis_products/by_type/spectrum_per_group_stack/<target>_<proposal>_src<NN>.png
  analysis_products/by_type/line_diagnostics/<line_name>/<target>_<proposal>_src<NN>_<line_name>.png

Each entry is a SYMLINK to the original; targets get the same source-of-truth
file. Re-running the script refreshes / re-creates the links (overwrites).
"""
import re
from pathlib import Path

ANALYSIS = Path("/orange/adamginsburg/salt/survey_2026/analysis_products")
BYTYPE = ANALYSIS / "by_type"

# pattern: source_<NN>_<line>_diagnostic.png   OR   <line>_diagnostic.png (legacy)
RE_DIAG_NEW = re.compile(r"^source_(\d{2})_(.+)_diagnostic\.png$")
RE_DIAG_OLD = re.compile(r"^(.+)_diagnostic\.png$")
RE_ALIGNED = re.compile(r"^aligned_by_(.+)\.png$")
RE_PANELS = re.compile(r"^spectrum_panels_src(\d{2})(?:_p(\d+))?\.png$")
RE_SALT = re.compile(r"^salt_stack\.png$")
RE_GROUPSTACK = re.compile(r"^spectrum_per_group_stack\.png$")


def make_link(src: Path, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    # number of '..' to climb from dest's directory back to ANALYSIS
    up = len(dest.parent.relative_to(ANALYSIS).parts)
    rel = Path(*([".."] * up)) / src.relative_to(ANALYSIS)
    dest.symlink_to(rel)


def main():
    n_overview = n_spec = n_kin = n_salt = n_grp = n_diag = 0
    for prop_dir in sorted(ANALYSIS.glob("*/2*")):
        if not prop_dir.is_dir():
            continue
        target = prop_dir.parent.name
        proposal = prop_dir.name

        # 1. overview_continuum.png
        ovc = prop_dir / "overview_continuum.png"
        if ovc.exists():
            make_link(ovc,
                       BYTYPE / "overview_continuum" /
                       f"{target}_{proposal}.png")
            n_overview += 1

        # 2. spectrum_panels_src<NN>[_pN].png
        for spec in prop_dir.glob("spectrum_panels_src*.png"):
            m = RE_PANELS.match(spec.name)
            if not m:
                continue
            sid = m.group(1); page = m.group(2)
            suffix = f"_p{page}" if page else ""
            make_link(spec,
                       BYTYPE / "spectrum_full" /
                       f"{target}_{proposal}_src{sid}{suffix}.png")
            n_spec += 1

        # 3. kinematic_stack/aligned_by_<guide>.png
        ksdir = prop_dir / "kinematic_stack"
        if ksdir.is_dir():
            for ks in ksdir.glob("aligned_by_*.png"):
                m = RE_ALIGNED.match(ks.name)
                if not m:
                    continue
                guide = m.group(1)
                make_link(ks,
                           BYTYPE / "kinematic_stack" /
                           f"{target}_{proposal}_{guide}.png")
                n_kin += 1

        # 4. salt_stack.png, spectrum_per_group_stack.png, and per-line
        #    *_diagnostic.png live in source_<NN>/
        for src_dir in sorted(prop_dir.glob("source_*")):
            m_src = re.match(r"^source_(\d{2})$", src_dir.name)
            if not m_src:
                continue
            sid = m_src.group(1)

            for sst in src_dir.glob("salt_stack.png"):
                make_link(sst,
                           BYTYPE / "salt_stack" /
                           f"{target}_{proposal}_src{sid}.png")
                n_salt += 1
            for gst in src_dir.glob("spectrum_per_group_stack.png"):
                make_link(gst,
                           BYTYPE / "spectrum_per_group_stack" /
                           f"{target}_{proposal}_src{sid}.png")
                n_grp += 1

            for diag in src_dir.glob("*_diagnostic.png"):
                m = RE_DIAG_NEW.match(diag.name)
                if m:
                    line = m.group(2)
                else:
                    m2 = RE_DIAG_OLD.match(diag.name)
                    if not m2:
                        continue
                    line = m2.group(1)
                line_clean = line.replace("/", "_")
                make_link(diag,
                           BYTYPE / "line_diagnostics" / line_clean /
                           f"{target}_{proposal}_src{sid}_{line_clean}.png")
                n_diag += 1

    print(f"linked overview_continuum: {n_overview}")
    print(f"linked spectrum_full:      {n_spec}")
    print(f"linked kinematic_stack:    {n_kin}")
    print(f"linked salt_stack:         {n_salt}")
    print(f"linked spectrum_per_group: {n_grp}")
    print(f"linked line_diagnostics:   {n_diag}")
    print(f"\nLine catalog (line_diagnostics subdirs):")
    for d in sorted((BYTYPE / "line_diagnostics").glob("*")):
        if d.is_dir():
            print(f"  {d.name}: {sum(1 for _ in d.iterdir())} sources")


if __name__ == "__main__":
    main()
