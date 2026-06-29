"""Bhula-style XCLASS runner adapted for our extracted source spectra.

Per (target_label):
  - locate spectra/<label>_spw*.fits and the corresponding <label>_band6.molfit
  - call task_myXCLASS.myXCLASS with the right freq grid, telescope size
    (= half-power beam), and v_LSR pulled from the target dir's metadata
  - write model spectrum + a data + model overlay PNG per spw to models/ and plots/

Activated env (Bhula): xclass_wrapper imports task_myXCLASS from
/orange/adamginsburg/software/XCLASS-Interface (via XCLASSRootDir env var).

Run via run_xclass.sh after `source ...` the xclass env.
"""
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits

HERE = Path(__file__).parent

# Pull XCLASS Python path from the env
XCLASSRootDir = os.environ.get("XCLASSRootDir",
                                "/orange/adamginsburg/software/XCLASS-Interface")
xclass_tasks = os.path.join(XCLASSRootDir, "build_tasks")
if xclass_tasks not in sys.path:
    sys.path.insert(0, xclass_tasks)

# Try a couple of import paths used by Bhula
try:
    from task_myXCLASS import myXCLASS
except ImportError as e:
    raise SystemExit(
        f"task_myXCLASS not importable from {xclass_tasks}. "
        f"Did you `source ...` the xclass env? Original: {e}"
    )


# Hand-curated v_LSR (km/s) per target label
VLSR_KMS = {
    "G326.6618+00.5207": -39.6,
    "NGC6334I": -7.0,
    "G015.0357-00.6795": 19.2,
    "G345.4938+01.4677": -12.6,
    "G345.5043+00.3480": -17.0,
    "MonR2-IRS3": 10.0,
    "Orion_SrcI": 5.0,
}

# Lines we want annotated on the overlay plots
ANNOTATE_GHZ = {
    "NaCl v=0 17-16": 217.97714756,
    "NaCl v=1 17-16": 219.6129951,
    "NaCl v=2 18-17": 230.7757972,
    "NaCl v=0 18-17": 234.25090148,
    "KCl v=0 28-27": 215.0083,
    "KCl v=0 30-29": 230.32056,
    "H2O 5(5,0)-6(4,3)": 232.6867,
    "H2O v2=1 3(1,3)-2(2,0)": 232.93657,
    "H2O v2=1 4(1,4)-3(2,1)": 234.1405,
    "H30alpha": 231.9009,
    "CH3CN 12-11 K=3": 220.7090,
    "CH3OH 8(0)-7(1) E": 220.0785,
    "OCS 18-17": 218.9034,
}


def label(spec_fits: Path) -> str:
    """Extract spw token from filename like G326_spw25.fits."""
    parts = spec_fits.stem.split("_")
    for t in parts:
        if t.startswith("spw"):
            return t
    return "spw??"


def run_xclass_one(spec_fits: Path, molfit_file: Path, vlsr: float,
                     model_dir: Path, plot_dir: Path):
    print(f"  -> {spec_fits.name}", flush=True)
    with fits.open(spec_fits) as hd:
        h = hd[0].header
        spec_data = np.asarray(hd[0].data).squeeze()
    nchan = spec_data.size
    nu0 = float(h["CRVAL1"])  # Hz
    dnu = float(h["CDELT1"])  # Hz
    freq_Hz = nu0 + np.arange(nchan) * dnu
    freq_MHz = freq_Hz / 1e6
    # XCLASS expects freq grid in MHz and FreqStep in MHz
    FreqMin = float(freq_MHz.min())
    FreqMax = float(freq_MHz.max())
    FreqStep = float(abs(dnu) / 1e6)
    bmaj = float(h.get("BMAJ", 0)) * 3600.0  # arcsec
    bmin = float(h.get("BMIN", 0)) * 3600.0
    tel_arcsec = 0.5 * (bmaj + bmin) if (bmaj > 0 and bmin > 0) else 0.5

    # Pass the spectrum file path to XCLASS so it loads + compares
    res = myXCLASS(
        FileName=str(spec_fits),
        FreqMin=FreqMin, FreqMax=FreqMax, FreqStep=FreqStep,
        TelescopeSize=tel_arcsec, Inter_Flag=True,
        t_back_flag=True, tBack=2.7, vLSR=vlsr,
        iso_flag=False,
        MolfitsFileBaseNames=[molfit_file.stem],
        IsoTableFileName="",
        CollisionFileName="",
        BackgroundFileName="",
        DustFileName="",
        verbose=False,
        Inter_Flag_Use=True,
        NumberProcessors=4,
    )
    # task_myXCLASS returns (modeldata, log, transitions, components, ...)
    modeldata = res[0]
    transitions = res[2] if len(res) > 2 else []
    # modeldata is 2D: [freq_MHz, T_K]
    mdl = np.asarray(modeldata)
    if mdl.ndim == 2 and mdl.shape[1] >= 2:
        mfreq_MHz = mdl[:, 0]; mT = mdl[:, 1]
    else:
        mfreq_MHz = freq_MHz; mT = np.full_like(freq_MHz, np.nan)
    spw = label(spec_fits)
    out_npz = model_dir / f"{spw}.model.npz"
    np.savez(out_npz, freq_MHz=mfreq_MHz, T_K=mT)

    # Diagnostic plot
    fig, ax = plt.subplots(figsize=(16, 4))
    ax.plot(freq_MHz / 1e3, spec_data, "k-", lw=0.6, label="data")
    ax.plot(mfreq_MHz / 1e3, mT, "C3-", lw=0.8, alpha=0.75,
             label="XCLASS LTE model")
    ax.set_xlabel("frequency (GHz)")
    ax.set_ylabel("T (K)")
    ax.set_title(f"{spec_fits.parent.parent.name}  {spw}  "
                  f"vLSR={vlsr:+.1f}")
    fmin = freq_MHz.min() / 1e3; fmax = freq_MHz.max() / 1e3
    ymax = float(np.nanpercentile(spec_data, 99)) * 1.15 + 1
    ax.set_ylim(float(np.nanpercentile(spec_data, 1)), ymax)
    for n, f in ANNOTATE_GHZ.items():
        if fmin <= f <= fmax:
            ax.axvline(f, color="C1", lw=0.4, alpha=0.6)
            ax.text(f, ymax * 0.97, n, rotation=90, ha="center", va="top",
                     fontsize=6, color="C1", alpha=0.85, clip_on=True)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(plot_dir / f"{spw}.png", dpi=120, bbox_inches="tight")
    plt.close(fig)


def main():
    if len(sys.argv) < 2:
        raise SystemExit("usage: xclass_runner.py <label>")
    label_name = sys.argv[1]
    tgt_dir = HERE / label_name
    spectra_dir = tgt_dir / "spectra"
    if not spectra_dir.is_dir():
        raise SystemExit(f"no spectra/ in {tgt_dir}")
    # Find molfit file
    molfits = list(tgt_dir.glob("*band6.molfit"))
    if not molfits:
        molfits = list(tgt_dir.glob("*.molfit"))
    if not molfits:
        raise SystemExit(f"no .molfit in {tgt_dir}")
    molfit = molfits[0]
    vlsr = VLSR_KMS.get(label_name, 0.0)
    print(f"label={label_name}, molfit={molfit.name}, v_LSR={vlsr}")

    model_dir = tgt_dir / "models"
    plot_dir = tgt_dir / "plots"
    model_dir.mkdir(exist_ok=True)
    plot_dir.mkdir(exist_ok=True)

    spec_fits_list = sorted(spectra_dir.glob(f"{label_name}_spw*.fits"))
    if not spec_fits_list:
        # fall back to any spw*.fits in spectra/
        spec_fits_list = sorted(spectra_dir.glob("*_spw*.fits"))
    if not spec_fits_list:
        raise SystemExit(f"no spw*.fits in {spectra_dir}")
    for f in spec_fits_list:
        try:
            run_xclass_one(f, molfit, vlsr, model_dir, plot_dir)
        except (RuntimeError, ValueError, KeyError) as e:
            print(f"    XCLASS error for {f.name}: {e}", flush=True)
    print(f"\nWrote models -> {model_dir}/")
    print(f"Wrote plots  -> {plot_dir}/")


if __name__ == "__main__":
    main()
