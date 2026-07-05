"""Run XCLASS LTE synthesis against each per-SPW source spectrum and overlay
the model on the data.

Modeled after /orange/adamginsburg/w51/hotcores/linecatalog/ethanbhula/xclass_wrapper.py
which calls task_myXCLASS.myXCLASS for the heavy lifting. Workflow:

  1. for each spectra/G326_<spw>.fits, read freq range
  2. invoke myXCLASS with G326_band6.molfit, T_back=2.7, vLSR=-39.6 km/s,
     TelescopeSize = mean(BMAJ, BMIN) in arcsec from the spectrum's FITS header
  3. write model spectrum to models/G326_<spw>.fits
  4. produce a 2-row diagnostic PNG: top = data + model overlay; bottom =
     residual; annotate the rest frequencies of NaCl 217.97/219.61/230.78/234.25,
     KCl rotational lines in band, and H2O 232.687 / 232.94 disk tracers.

XCLASS lives under XCLASSRootDir; set XCLASSRootDir env var or rely on the
default /orange/adamginsburg/software/XCLASS-Interface .

This is a runner stub: actually invoking XCLASS requires the XCLASS Python
environment activated (see Bhula's README). Run inside that env.
"""
import os
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

HERE = Path(__file__).parent
XCLASSRootDir = os.environ.get("XCLASSRootDir",
                                 "/orange/adamginsburg/software/XCLASS-Interface")
sys.path.insert(0, os.path.join(XCLASSRootDir, "build_tasks"))

MOLFIT = HERE / "G326_band6.molfit"
SPECTRA = HERE / "spectra"
MODELS = HERE / "models"
PLOTS = HERE / "plots"
MODELS.mkdir(exist_ok=True)
PLOTS.mkdir(exist_ok=True)

V_LSR = -39.6
T_BACK = 2.7
NUM_PROCESSORS = 8

# Lines of interest for line-ID annotation (rest GHz)
ANNOTATE = {
    "NaCl v=0 J=17-16": 217.97714756,
    "NaCl v=1 J=17-16": 219.6129951,
    "NaCl v=2 J=18-17": 230.7757972,
    "NaCl v=0 J=18-17": 234.25090148,
    "KCl v=0 J=29-28": 222.0,  # placeholder; lookup exact
    "H2O 5(5,0)-6(4,3) 232.69": 232.6867,
    "H2O v2=1 3(1,3)-2(2,0) 232.94": 232.93657,
    "H30alpha": 231.90,
}


def run_xclass(spw_fits: Path):
    try:
        from task_myXCLASS import myXCLASS
    except ImportError:
        print(f"XCLASS not importable (XCLASSRootDir={XCLASSRootDir}). "
               f"Source the XCLASS environment and rerun.")
        return None
    with fits.open(spw_fits) as hd:
        h = hd[0].header
        spec = hd[0].data.squeeze()
    freq_lo = float(h["CRVAL1"]) / 1e6  # MHz
    dnu = float(h["CDELT1"]) / 1e6  # MHz
    nchan = spec.size
    freq_MHz = freq_lo + np.arange(nchan) * dnu
    TelescopeSize = 0.5 * (float(h.get("BMAJ", 0)) + float(h.get("BMIN", 0))) * 3600.0

    out = myXCLASS(
        FreqMin=float(freq_MHz.min()),
        FreqMax=float(freq_MHz.max()),
        FreqStep=float(abs(dnu)),
        TelescopeSize=TelescopeSize,
        Inter_Flag=True,
        t_back_flag=True,
        tBack=T_BACK,
        vLSR=V_LSR,
        MolfitsFileName=str(MOLFIT),
        IsoTableFileName="",
        NumberProcessors=NUM_PROCESSORS,
    )
    return out


def main():
    for spec_fits in sorted(SPECTRA.glob("G326_spw*.fits")):
        print(f"==> {spec_fits.name}")
        result = run_xclass(spec_fits)
        if result is None:
            continue
        # XCLASS returns (modeldata, log, transitions, components, ...)
        modeldata = result[0]
        np.savetxt(MODELS / spec_fits.with_suffix(".dat").name,
                    modeldata)
        print(f"  wrote {MODELS / spec_fits.stem}.dat")


if __name__ == "__main__":
    main()
