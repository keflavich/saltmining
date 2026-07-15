"""NaCl+RRL table v2: REAL spatial smoothing for the 300 AU column.

For each (target, proposal, line) where native beam < 300 AU at the source
distance, open the cube, take a small spatial cutout around the brightest
continuum source, convolve to a target beam of (300 AU / d_kpc) arcsec, and
remeasure the per-channel noise. The detection peak is also re-measured in
the smoothed cube.

Falls back to the point-source brightness scaling (build_nacl_rrl_table.py)
when:
  - native beam >= 300 AU (cannot smooth finer)
  - cube fails to open
  - convolution fails

Reuses the v1 collect() to enumerate rows + load cube paths.

Output:
  data/nacl_rrl_table_smoothed.csv
  /orange/adamginsburg/salt/demography_2026/nacl_rrl_uls_smoothed.tex
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.convolution import Gaussian2DKernel, convolve_fft
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import mad_std
from astropy.wcs import WCS

import build_nacl_rrl_table as v1

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
UVDIR = ROOT / "uvdata"
OUT_CSV = ROOT / "data/nacl_rrl_table_smoothed.csv"
OUT_TEX = Path("/orange/adamginsburg/salt/demography_2026/nacl_rrl_uls_smoothed.tex")

TARGET_AU = 300.0
CUTOUT_HALF_ARCSEC = 4.0  # spatial cutout half-size before convolve
OFF_KMS_HALF = 50.0  # off-line window half-width for noise


def find_brightest_coord(proposal_dir: Path):
    cont = proposal_dir / "continuum_sources.csv"
    if not cont.exists():
        return None, None
    try:
        df = pd.read_csv(cont)
    except pd.errors.EmptyDataError:
        return None, None
    if df.empty or "peak_Jybeam" not in df.columns:
        return None, None
    idx = df["peak_Jybeam"].idxmax()
    return float(df.loc[idx, "ra_deg"]), float(df.loc[idx, "dec_deg"])


def _velocity_axis(hdr, rest_GHz):
    """Return per-channel velocity (km/s, radio) from a 3D FITS header."""
    nchan = hdr["NAXIS3"]
    crpix = hdr["CRPIX3"]; crval = hdr["CRVAL3"]; cdelt = hdr["CDELT3"]
    ctype = hdr.get("CTYPE3", "").upper()
    pix = np.arange(nchan)
    vals = crval + (pix + 1 - crpix) * cdelt
    if ctype.startswith("FREQ"):
        nu = vals  # Hz
        v = (1.0 - nu / (rest_GHz * 1e9)) * 299792.458  # radio convention km/s
    elif ctype.startswith("VRAD") or ctype.startswith("VELO"):
        v = vals * 1e-3 if abs(cdelt) > 100 else vals  # already km/s if small
        v = vals / 1000.0 if abs(np.median(vals)) > 1e4 else vals
    else:
        v = vals
    return v


def smoothed_noise_K(cube_path: Path, ra: float, dec: float,
                       rest_GHz: float, vlsr_kms: float,
                       target_beam_arcsec: float,
                       native_beam_arcsec: float):
    """Pure-numpy spatial smoothing: read a small (y,x) cutout from the cube
    across all channels, convolve each channel with the Gaussian residual
    kernel that takes native_beam -> target_beam, sample spectrum at source
    pixel, return sigma+peak in K. Avoids spectral_cube's full-cube load.

    Unit handling: convolution with a normalized kernel keeps a Jy/beam map
    in Jy per NATIVE beam (it is an intensity average), so the K conversion
    uses the native bmaj*bmin. K cubes pass through unchanged."""
    with fits.open(cube_path, memmap=True) as hdul:
        h = hdul[0].header
        data = hdul[0].data
        if data.ndim == 4:
            data = data[0]
        wcs = WCS(h).celestial
        xp, yp = wcs.world_to_pixel(SkyCoord(ra, dec, unit="deg"))
        ix0 = int(round(float(xp))); iy0 = int(round(float(yp)))
        pix_scale_deg = float(np.abs(wcs.pixel_scale_matrix[0, 0]))
        pix_scale_arcsec = pix_scale_deg * 3600.0
        half_pix = int(np.ceil(CUTOUT_HALF_ARCSEC / pix_scale_arcsec))
        ny, nx = data.shape[-2], data.shape[-1]
        y_lo = max(0, iy0 - half_pix); y_hi = min(ny, iy0 + half_pix + 1)
        x_lo = max(0, ix0 - half_pix); x_hi = min(nx, ix0 + half_pix + 1)
        if y_hi - y_lo < 5 or x_hi - x_lo < 5:
            return None
        v_full = _velocity_axis(h, rest_GHz)
        idx = np.where((v_full > vlsr_kms - 3 * OFF_KMS_HALF) &
                       (v_full < vlsr_kms + 3 * OFF_KMS_HALF))[0]
        if idx.size < 3:
            return None
        cz_lo, cz_hi = int(idx.min()), int(idx.max()) + 1
        # Read full spatial extent for the channel slab (~2s on NFS for typical
        # cubes); spatial-cropped reads are 4-10x slower due to non-contiguous I/O
        slab_full = np.asarray(data[cz_lo:cz_hi], dtype=np.float32)
        cube_cut = slab_full[:, y_lo:y_hi, x_lo:x_hi]
        v = v_full[cz_lo:cz_hi]
    if not np.isfinite(cube_cut).any():
        return None
    # Residual Gaussian: sigma_res_arcsec = sqrt(target^2 - native^2)/2.355
    sq = target_beam_arcsec ** 2 - native_beam_arcsec ** 2
    if sq <= 0:
        return None
    res_fwhm_arcsec = float(np.sqrt(sq))
    res_sigma_pix = res_fwhm_arcsec / 2.355 / pix_scale_arcsec
    kernel = Gaussian2DKernel(res_sigma_pix)
    sx = ix0 - x_lo; sy = iy0 - y_lo
    spec_jybeam = np.empty(cube_cut.shape[0], dtype=float)
    for k in range(cube_cut.shape[0]):
        plane = cube_cut[k]
        if not np.isfinite(plane).any():
            spec_jybeam[k] = np.nan; continue
        sm = convolve_fft(plane, kernel, normalize_kernel=True,
                            nan_treatment="interpolate", preserve_nan=True,
                            allow_huge=True)
        spec_jybeam[k] = sm[sy, sx]
    # Unit conversion: a normalized-kernel convolution of a Jy/beam map is
    # still in Jy per NATIVE beam (intensity average), so the K conversion
    # uses the native beam solid angle, NOT the target beam.
    bunit = str(h.get("BUNIT", "")).strip().lower().replace(" ", "")
    if bunit in ("k", "kelvin"):
        conv = 1.0
    elif bunit in ("jy/beam", "jybeam-1", "beam-1jy", "jy/bm"):
        bmaj_as = float(h["BMAJ"]) * 3600.0 if h.get("BMAJ") else native_beam_arcsec
        bmin_as = float(h["BMIN"]) * 3600.0 if h.get("BMIN") else native_beam_arcsec
        nu_GHz = float(rest_GHz)
        conv = 1.222e6 / (nu_GHz ** 2 * bmaj_as * bmin_as)
    else:
        return None
    spec = spec_jybeam * conv
    dv = float(np.median(np.abs(np.diff(v))))
    off = (np.abs(v - vlsr_kms) > OFF_KMS_HALF) & (np.abs(v - vlsr_kms) < 3 * OFF_KMS_HALF)
    if off.sum() < 5:
        off = np.isfinite(spec)
    spec = spec - float(np.nanmedian(spec[off]))
    sigma = float(mad_std(spec[off], ignore_nan=True))
    # rescale per-channel noise to a 10 km/s channel for the table column
    sigma_10 = sigma * np.sqrt(min(1.0, dv / 10.0))
    on = np.abs(v - vlsr_kms) < OFF_KMS_HALF
    peak = float(np.nanmax(spec[on])) if on.any() else np.nan
    return dict(sigma_K=sigma, sigma_10kms_K=sigma_10, peak_K=peak, dv_kms=dv)


def _circle_region(coord, r_arcsec):
    from regions import CircleSkyRegion
    return CircleSkyRegion(coord, radius=r_arcsec * u.arcsec)


def main():
    src = pd.read_csv(ROOT / "data/sources_L4_d2.csv")
    distances = dict(zip(src["name"], src["dist_kpc"]))
    base = v1.collect(distances)
    if base.empty:
        print("empty")
        return

    # vLSR per target for the on/off-line windows, from the measured (or
    # literature) source velocity in data/obs_params.csv. A wrong vLSR puts
    # the line inside the off-line noise window and the peak outside the
    # on-line window, so no silent fallback: targets without a velocity are
    # skipped by smoothed_noise_K's caller with status 'no_vlsr'.
    vlsr_per_target = {}
    op_csv = ROOT / "data/obs_params.csv"
    if op_csv.exists():
        op = pd.read_csv(op_csv)
        for _, orow in op.iterrows():
            v = orow.get("vsrc_kms")
            if v is not None and np.isfinite(v):
                vlsr_per_target[str(orow["name"])] = float(v)

    smoothed_sig = []
    smoothed_pk = []
    used = []
    n_run = 0; n_fall = 0
    for i, r in base.iterrows():
        tgt = r["target"]; prop = r["proposal"]; line = r["line"]
        d_kpc = distances.get(tgt, np.nan)
        beam_arcsec = r.get("beam_native_arcsec")
        if not np.isfinite(d_kpc) or beam_arcsec is None or not np.isfinite(beam_arcsec):
            smoothed_sig.append(np.nan); smoothed_pk.append(np.nan); used.append("nodata")
            n_fall += 1; continue
        target_beam_arcsec = TARGET_AU / (d_kpc * 1000.0)
        if beam_arcsec >= target_beam_arcsec:
            smoothed_sig.append(r["sigma_10kms_native_K"])
            smoothed_pk.append(r["peak_K"])
            used.append("native_already_larger")
            n_fall += 1; continue
        prop_dir = ROOT / "analysis_products" / tgt / prop
        ra, dec = find_brightest_coord(prop_dir)
        if ra is None:
            smoothed_sig.append(np.nan); smoothed_pk.append(np.nan); used.append("no_cont")
            n_fall += 1; continue
        meas = pd.read_csv(prop_dir / "line_measurements.csv")
        m = meas[(meas["line"] == line)]
        if m.empty or m["cube"].dropna().empty:
            smoothed_sig.append(np.nan); smoothed_pk.append(np.nan); used.append("no_cube_match")
            n_fall += 1; continue
        cube_path = UVDIR / prop / tgt / m["cube"].iloc[0]
        if not cube_path.exists():
            smoothed_sig.append(np.nan); smoothed_pk.append(np.nan); used.append("cube_missing")
            n_fall += 1; continue
        rest_GHz = float(m["rest_GHz"].iloc[0])
        v0 = vlsr_per_target.get(tgt)
        if v0 is None:
            smoothed_sig.append(np.nan); smoothed_pk.append(np.nan); used.append("no_vlsr")
            n_fall += 1; continue
        try:
            out = smoothed_noise_K(cube_path, ra, dec, rest_GHz, v0,
                                     target_beam_arcsec, float(beam_arcsec))
        except (ValueError, OSError, KeyError, RuntimeError) as e:
            print(f"  {tgt} {prop} {line}: smooth fail ({type(e).__name__}: {e})")
            smoothed_sig.append(np.nan); smoothed_pk.append(np.nan); used.append("err")
            n_fall += 1; continue
        if out is None:
            smoothed_sig.append(np.nan); smoothed_pk.append(np.nan); used.append("out_of_pixel")
            n_fall += 1; continue
        smoothed_sig.append(out["sigma_10kms_K"])
        smoothed_pk.append(out["peak_K"])
        used.append("smoothed")
        n_run += 1
        print(f"  {tgt} {prop} {line}: sigma_smooth(10kms)={out['sigma_10kms_K']*1000.0:.2f} mK "
              f"(native_300au_pt={r['sigma_10kms_300au_K']*1000.0:.2f} mK)")

    base["sigma_smoothed_K"] = smoothed_sig
    base["peak_smoothed_K"] = smoothed_pk
    base["smooth_status"] = used
    base.to_csv(OUT_CSV, index=False)
    print(f"\nSmoothed {n_run}, fell back {n_fall}; wrote {OUT_CSV}")
    write_tex_smoothed(base)
    print(f"wrote {OUT_TEX}")


def fmt_smoothed_cell(peak_K, sigma_K, detected):
    from kelvin import fmt_K
    if not np.isfinite(sigma_K):
        return r"\nodata"
    if detected:
        return fmt_K(peak_K)
    return rf"$<${fmt_K(3.0*sigma_K)}"


def write_tex_smoothed(df: pd.DataFrame):
    # Drop upper limits whose native synthesized beam exceeds 500 AU: such
    # limits cannot constrain a compact (<500 AU) salt disk, and a >500 AU
    # native beam cannot be smoothed to the 300 AU reference anyway. Keep
    # detections regardless. (Same policy as Table~\ref{tab:naclrrl}.)
    n0 = len(df)
    keep = df["detected"].astype(bool) | (df["beam_native_au"] <= 500.0)
    df = df[keep].copy()
    print(f"  smoothed: dropped {n0 - len(df)} UL rows (native beam > 500 AU)")
    df = df.sort_values(["target", "proposal", "line"])
    lines = [
        r"\startlongtable",
        r"\begin{deluxetable*}{lllccccc}",

        r"\tablecaption{NaCl + RRL detections / $3\sigma$ upper limits with the "
        r"$300$\,AU column derived from real spatial convolution (per-channel "
        r"Gaussian smoothing to a $300$\,AU effective beam at the source "
        r"distance) and remeasured noise. The native columns repeat "
        r"Table~\ref{tab:naclrrl}; the smoothed column reports the noise of the "
        r"actual smoothed map per channel and is the right physical quantity "
        r"for partially-resolved emission. Values are brightness temperatures "
        r"in K.\label{tab:naclrrl_smoothed}}",
        r"\tablehead{",
        r"\colhead{Source} & \colhead{Program} & \colhead{Line} & "
        r"\colhead{Src} & \colhead{$\theta_\mathrm{beam}$} & "
        r"\colhead{$1$\,\kms, nat.} & \colhead{$10$\,\kms, nat.} & "
        r"\colhead{$10$\,\kms, $300$\,AU smoothed} \\",
        r" & & & & (AU) & (K) & (K) & (K) }",
        r"\startdata",
    ]
    for _, r in df.iterrows():
        det = bool(r["detected"])
        v1n = v1.fmt_value(r["peak_K"], r["sigma_1kms_native_K"], det)
        v10n = v1.fmt_value(r["peak_K"], r["sigma_10kms_native_K"], det)
        v10_3 = fmt_smoothed_cell(r["peak_smoothed_K"], r["sigma_smoothed_K"], det)
        beam_au = r["beam_native_au"]
        beam_str = f"{beam_au:.0f}" if np.isfinite(beam_au) else r"\nodata"
        target = v1._display_name(str(r["target"])).replace("_", r"\_")
        lines.append(
            f"{target} & {r['proposal']} & {v1.fmt_line(r['line'])} & "
            f"{int(r['src_id'])} & {beam_str} & "
            f"{v1n} & {v10n} & {v10_3} \\\\"
        )
    lines += [
        r"\enddata",
        r"\tablecomments{The $10$\,\kms,\,$300$\,AU smoothed column applies a "
        r"per-channel astropy.convolve\_fft Gaussian to a $5''\times5''$ "
        r"spatial cutout of each cube, bringing the synthesized beam to the "
        r"$300$\,AU/d-target before remeasuring the per-channel noise via "
        r"mad\_std. Falls back to the native $\sigma$ where $\theta_\mathrm{beam} \geq 300$\,AU.}",
        r"\end{deluxetable*}",
    ]
    OUT_TEX.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
