"""Build per-source diagnostic moment-map panels for the 'first reported
here' sources in Table 4.

For each source we render a single multi-panel PNG combining the brightest
NaCl, KCl, H2O 232 GHz, and SiS/SO mom0 maps (up to 4 panels). The PNGs
are saved next to the source's analysis_products dir under a `figures/`
subdir for inclusion in main.tex.

Targets are determined from data_summary's footnote system (any row with
* footnote, i.e., new this work).
"""
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ImageNormalize, AsymmetricPercentileInterval

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
OUT_FIG = Path("/orange/adamginsburg/salt/demography_2026/figures")
OUT_FIG.mkdir(parents=True, exist_ok=True)

# Map display label (Table-4 row) -> directory name on disk.
NEW_DETECTION_TARGETS = [
    ("IRAS 15412-5359",  "G326.6618+00.5207"),
    ("IRAS 18174-1612",  "G015.0357-00.6795"),
    ("I17016-4124",      "G345.5043+00.3480"),
    ("IRAS 17233-3606",  "IRAS17233-3606"),
]

# Per source, choose up to 4 mom0 maps to render. Priority order.
LINE_PRIORITIES = [
    ("NaCl", r"NaCl"),
    ("KCl",  r"KCl"),
    ("H2O 232", r"H2O.*232"),
    ("SiS",  r"SiS"),
    ("SO",   r"SO_"),
    ("H2O", r"H2O_"),
]


def best_proposal_dir(target_dir):
    """Pick the proposal subdir with the most line_measurements."""
    best = None
    best_n = -1
    for p in sorted(target_dir.glob("2*")):
        meas = p / "line_measurements.csv"
        if not meas.exists():
            continue
        try:
            n = sum(1 for _ in meas.open()) - 1
        except OSError:
            n = 0
        if n > best_n:
            best_n = n
            best = p
    return best


def brightest_source(prop_dir):
    cont = prop_dir / "continuum_sources.csv"
    if not cont.exists():
        return None
    import pandas as pd
    try:
        df = pd.read_csv(cont)
    except pd.errors.EmptyDataError:
        return None
    if df.empty or "peak_Jybeam" not in df.columns:
        return None
    return int(df.loc[df["peak_Jybeam"].idxmax(), "id"])


def find_mom0(prop_dir, src_id, line_label, pattern):
    """Return the highest-priority mom0 fits matching pattern at src_id."""
    import re
    src_dir = prop_dir / f"source_{src_id:02d}"
    if not src_dir.is_dir():
        return None
    rx = re.compile(pattern, re.IGNORECASE)
    candidates = [p for p in src_dir.glob("*_mom0.fits") if rx.search(p.name)]
    # Sort: prefer v=0 (ground state, less risk of confusion) for NaCl/KCl,
    # else first.
    def key(p):
        n = p.name
        v0 = 0 if "_v0_" in n else 1
        return (v0, n)
    candidates.sort(key=key)
    return candidates[0] if candidates else None


def render_panel(ax, fits_path, title):
    with fits.open(fits_path) as hdul:
        data = hdul[0].data
        hdr = hdul[0].header
    if data is None:
        ax.set_axis_off()
        ax.set_title(f"{title} (no data)")
        return
    while data.ndim > 2:
        data = data[0]
    wcs = WCS(hdr).celestial
    norm = ImageNormalize(data, interval=AsymmetricPercentileInterval(2, 99.5))
    im = ax.imshow(data, origin="lower", cmap="inferno", norm=norm)
    ax.set_title(title, fontsize=10)
    # Strip ticks for compactness
    ax.tick_params(axis="both", which="both", labelsize=7)


def build_figure(display, target):
    tdir = ANALYSIS / target
    if not tdir.is_dir():
        for alias in (target.replace("-", "_"), target.replace("_", "-")):
            if (ANALYSIS / alias).is_dir():
                tdir = ANALYSIS / alias
                break
    if not tdir.is_dir():
        print(f"  no analysis dir for {target}")
        return None
    pd_ = best_proposal_dir(tdir)
    if pd_ is None:
        print(f"  no proposal dir under {tdir}")
        return None
    bid = brightest_source(pd_)
    if bid is None:
        print(f"  no brightest source in {pd_}")
        return None
    panels = []
    for label, pat in LINE_PRIORITIES:
        m = find_mom0(pd_, bid, label, pat)
        if m is not None:
            panels.append((label, m))
        if len(panels) == 4:
            break
    if not panels:
        print(f"  no mom0 maps for {target} src_{bid:02d}")
        return None
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4), squeeze=False)
    for ax, (label, fp) in zip(axes[0], panels):
        render_panel(ax, fp, f"{label}: {fp.stem.replace('source_'+f'{bid:02d}_','')}")
    fig.suptitle(f"{display} (src_{bid:02d}, {pd_.name})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    safe = target.replace("/", "_").replace(" ", "_")
    outp = OUT_FIG / f"newdet_{safe}_mom0.png"
    fig.savefig(outp, dpi=120)
    plt.close(fig)
    print(f"  wrote {outp.name}  ({n} panels)")
    return outp


def main():
    for display, target in NEW_DETECTION_TARGETS:
        print(f"== {display}  ({target})")
        build_figure(display, target)


if __name__ == "__main__":
    main()
