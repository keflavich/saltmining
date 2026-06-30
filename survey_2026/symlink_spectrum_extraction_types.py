"""Add suffix-disambiguated symlinks for existing spectrum files so it's
clear at a glance whether you're looking at a peak-pixel, aperture-mean,
or kinematic-stack extraction.

Mapping:
  spectra/lineid_full.{png,pdf}             -> spectra/peakpixel_lineid_full.{png,pdf}
  source_NN/spectrum_panels_src*.{png,pdf}  -> source_NN/aperture_spectrum_panels_src*.{png,pdf}
  kinematic_stack/aligned_by_*.{png,npz}    -> already named explicitly; nothing to do

Doesn't delete originals -- the symlinks are pure additions so existing
references still work."""
from pathlib import Path

ANALYSIS = Path("/orange/adamginsburg/salt/survey_2026/analysis_products")


def main():
    n_added = 0
    for tdir in sorted(ANALYSIS.iterdir()):
        if not tdir.is_dir() or tdir.name.startswith("_"):
            continue
        # spectra/ at target level
        sdir = tdir / "spectra"
        if sdir.is_dir():
            for ext in ("png", "pdf"):
                src = sdir / f"lineid_full.{ext}"
                dst = sdir / f"peakpixel_lineid_full.{ext}"
                if src.exists() and not dst.exists():
                    try:
                        dst.symlink_to(src.name)
                        n_added += 1
                    except OSError:
                        pass
        # spectrum_panels under each proposal's source_NN/
        for pdir in tdir.glob("2*"):
            for sdir2 in pdir.glob("source_*"):
                if not sdir2.is_dir():
                    continue
                for f in list(sdir2.glob("spectrum_panels_src*.png")) + \
                          list(sdir2.glob("spectrum_panels_src*.pdf")):
                    dst = f.parent / f"aperture_{f.name}"
                    if not dst.exists():
                        try:
                            dst.symlink_to(f.name)
                            n_added += 1
                        except OSError:
                            pass
    print(f"added {n_added} symlinks")


if __name__ == "__main__":
    main()
