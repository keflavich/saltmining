"""Extract source-pixel spectra from all G326.6618+00.5207 cubes for XCLASS.

XCLASS expects ASCII / FITS spectrum files: one per SPW with the brightness
temperature at the YSO continuum peak as a function of frequency. We write
one FITS per SPW with axes (1, NCHAN) of K (Rayleigh-Jeans converted from
Jy/beam using the cube's beam) plus a 1D table for plotting.

Source = brightest continuum peak in analysis_products/.../continuum_sources.csv.

Each cube gets written as:
    spectra/G326_<spwid>.fits   (FITS, 1 axis spectral, BUNIT=K)
"""
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
TARGET = "G326.6618+00.5207"
PROP = "2022.1.01344.S"
UVDIR = ROOT / "uvdata" / PROP / TARGET
OUTDIR = Path(__file__).parent / "spectra"
OUTDIR.mkdir(exist_ok=True)


def jy_per_beam_to_K(nu_GHz, bmaj_arcsec, bmin_arcsec):
    return 1.222e6 / (nu_GHz ** 2 * bmaj_arcsec * bmin_arcsec)


def main():
    cont = pd.read_csv(ROOT / "analysis_products" / TARGET / PROP / "continuum_sources.csv")
    bid = int(cont.loc[cont["peak_Jybeam"].idxmax(), "id"])
    ra = float(cont.loc[bid - 1, "ra_deg"])
    dec = float(cont.loc[bid - 1, "dec_deg"])
    coord = SkyCoord(ra, dec, unit="deg")
    print(f"brightest source #{bid} at {ra:.6f}, {dec:.6f}")

    cubes = sorted(UVDIR.glob("*sci.spw*.cube.I.pbcor.fits"))
    print(f"found {len(cubes)} cubes")

    summary = []
    for cube_path in cubes:
        with fits.open(cube_path, memmap=True) as hdul:
            h = hdul[0].header
            data = hdul[0].data
        if data.ndim == 4:
            data = data[0]
        wcs = WCS(h).celestial
        xp, yp = wcs.world_to_pixel(coord)
        ix = int(round(float(xp))); iy = int(round(float(yp)))
        if not (0 <= ix < data.shape[2] and 0 <= iy < data.shape[1]):
            print(f"  {cube_path.name}: source off-FOV, skipping")
            continue
        spec_jybeam = data[:, iy, ix]
        nu0 = float(h.get("CRVAL3"))  # Hz
        dnu = float(h.get("CDELT3"))  # Hz
        crpix3 = float(h.get("CRPIX3"))
        nchan = data.shape[0]
        freq_Hz = nu0 + (np.arange(nchan) + 1 - crpix3) * dnu
        bmaj = float(h.get("BMAJ", 0)) * 3600
        bmin = float(h.get("BMIN", 0)) * 3600
        nu_center_GHz = float(np.median(freq_Hz)) / 1e9
        Kfac = jy_per_beam_to_K(nu_center_GHz, bmaj, bmin)
        spec_K = spec_jybeam * Kfac
        # XCLASS reads ASCII tables: freq_MHz  T_K  d_T_K
        freq_MHz = freq_Hz / 1e6
        spw = ""
        for tok in cube_path.stem.split("."):
            if tok.startswith("spw"):
                spw = tok; break
        # save FITS 2D (1, NCHAN) per XCLASS conventions
        out_fits = OUTDIR / f"G326_{spw}.fits"
        hdr = fits.Header()
        hdr["BUNIT"] = "K"
        hdr["CTYPE1"] = "FREQ"; hdr["CRPIX1"] = 1
        hdr["CRVAL1"] = float(freq_Hz[0]); hdr["CDELT1"] = float(dnu)
        hdr["CUNIT1"] = "Hz"; hdr["RESTFRQ"] = nu0
        hdr["BMAJ"] = float(h.get("BMAJ", 0)); hdr["BMIN"] = float(h.get("BMIN", 0))
        hdr["BPA"] = float(h.get("BPA", 0))
        hdr["OBJECT"] = TARGET
        hdr["SRC_RA"] = ra; hdr["SRC_DEC"] = dec
        hdr["SRC_ID"] = bid
        fits.PrimaryHDU(spec_K[np.newaxis, :].astype(np.float32), hdr).writeto(
            out_fits, overwrite=True)
        # ASCII table for XCLASS (3-column: freq_MHz, T_K, sigma_K)
        sigma_K = float(np.nanstd(spec_K[(spec_K < np.nanpercentile(spec_K, 90)) &
                                          (spec_K > np.nanpercentile(spec_K, 10))]))
        ascii_out = OUTDIR / f"G326_{spw}.dat"
        with open(ascii_out, "w") as fh:
            for k in range(nchan):
                fh.write(f"{freq_MHz[k]:.6f}  {spec_K[k]:.4f}  {sigma_K:.4f}\n")
        summary.append({
            "spw": spw, "fits": out_fits.name, "ascii": ascii_out.name,
            "nchan": nchan,
            "freq_lo_GHz": float(freq_Hz.min() / 1e9),
            "freq_hi_GHz": float(freq_Hz.max() / 1e9),
            "dnu_kHz": dnu / 1e3,
            "T_max_K": float(np.nanmax(spec_K)),
            "T_med_K": float(np.nanmedian(spec_K)),
            "sigma_K": sigma_K,
        })
        print(f"  {spw}: nu={summary[-1]['freq_lo_GHz']:.3f}-{summary[-1]['freq_hi_GHz']:.3f} GHz "
              f"T_max={summary[-1]['T_max_K']:.1f} K Tmed={summary[-1]['T_med_K']:.1f} K sigma={sigma_K:.2f} K")
    pd.DataFrame(summary).to_csv(OUTDIR.parent / "spw_summary.csv", index=False)
    print(f"\nwrote {OUTDIR}/ ({len(summary)} spectra)")


if __name__ == "__main__":
    main()
