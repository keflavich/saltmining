"""Extract source-pixel spectra from all cubes of a (target, proposal) for
XCLASS line ID. Generalized version of the G326-specific script.

For each cube:
  - source pixel = WCS-projected coordinate of the brightest mm continuum
    source (peak_Jybeam max in continuum_sources.csv)
  - extract spectrum at that pixel
  - convert Jy/beam -> K via per-cube beam + median frequency
  - write FITS (1, NCHAN, BUNIT=K) and 3-col ASCII (freq_MHz, T_K, sigma_K)

Usage:
    python extract_source_spectra.py --target G326.6618+00.5207 \
        --proposal 2022.1.01344.S --label G326

Output goes to lineid/<label>/spectra/.
"""
import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")


def jy_per_beam_to_K(nu_GHz, bmaj_arcsec, bmin_arcsec):
    return 1.222e6 / (nu_GHz ** 2 * bmaj_arcsec * bmin_arcsec)


def extract(target, proposal, label, source_id=None):
    uvdir = ROOT / "uvdata" / proposal / target
    if not uvdir.is_dir():
        raise SystemExit(f"no uvdir: {uvdir}")
    cont = pd.read_csv(ROOT / "analysis_products" / target / proposal /
                        "continuum_sources.csv")
    if source_id is None:
        bid = int(cont.loc[cont["peak_Jybeam"].idxmax(), "id"])
    else:
        bid = source_id
    ra = float(cont.loc[bid - 1, "ra_deg"])
    dec = float(cont.loc[bid - 1, "dec_deg"])
    coord = SkyCoord(ra, dec, unit="deg")
    print(f"target={target} proposal={proposal} src #{bid} at "
          f"{ra:.6f},{dec:.6f}")

    outdir = Path(__file__).parent / label / "spectra"
    outdir.mkdir(parents=True, exist_ok=True)

    cubes = sorted(uvdir.glob("*sci*cube*pbcor.fits"))
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
        nu0 = float(h.get("CRVAL3"))
        dnu = float(h.get("CDELT3"))
        crpix3 = float(h.get("CRPIX3"))
        nchan = data.shape[0]
        freq_Hz = nu0 + (np.arange(nchan) + 1 - crpix3) * dnu
        bmaj = float(h.get("BMAJ", 0)) * 3600
        bmin = float(h.get("BMIN", 0)) * 3600
        nu_center_GHz = float(np.median(freq_Hz)) / 1e9
        K = jy_per_beam_to_K(nu_center_GHz, bmaj, bmin)
        spec_K = spec_jybeam * K
        freq_MHz = freq_Hz / 1e6
        spw = next((t for t in cube_path.stem.split(".") if t.startswith("spw")),
                    "spw??")
        out_fits = outdir / f"{label}_{spw}.fits"
        hdr = fits.Header()
        hdr["BUNIT"] = "K"
        hdr["CTYPE1"] = "FREQ"; hdr["CRPIX1"] = 1
        hdr["CRVAL1"] = float(freq_Hz[0]); hdr["CDELT1"] = float(dnu)
        hdr["CUNIT1"] = "Hz"; hdr["RESTFRQ"] = nu0
        hdr["BMAJ"] = float(h.get("BMAJ", 0)); hdr["BMIN"] = float(h.get("BMIN", 0))
        hdr["BPA"] = float(h.get("BPA", 0))
        hdr["OBJECT"] = target
        hdr["SRC_RA"] = ra; hdr["SRC_DEC"] = dec; hdr["SRC_ID"] = bid
        fits.PrimaryHDU(spec_K[np.newaxis, :].astype(np.float32),
                          hdr).writeto(out_fits, overwrite=True)
        T_low = np.nanpercentile(spec_K, 10)
        T_hi = np.nanpercentile(spec_K, 90)
        sigma_K = float(np.nanstd(spec_K[(spec_K > T_low) & (spec_K < T_hi)]))
        ascii_out = outdir / f"{label}_{spw}.dat"
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
        print(f"  {spw}: nu={summary[-1]['freq_lo_GHz']:.3f}-"
              f"{summary[-1]['freq_hi_GHz']:.3f} GHz "
              f"T_max={summary[-1]['T_max_K']:.1f} K "
              f"T_med={summary[-1]['T_med_K']:.1f} K "
              f"sigma={sigma_K:.2f} K")
    pd.DataFrame(summary).to_csv(outdir.parent / "spw_summary.csv", index=False)
    print(f"wrote {outdir} ({len(summary)} spectra)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--proposal", required=True)
    ap.add_argument("--label", required=True,
                     help="output subdir under lineid/")
    ap.add_argument("--source-id", type=int, default=None)
    args = ap.parse_args()
    extract(args.target, args.proposal, args.label, args.source_id)


if __name__ == "__main__":
    main()
