"""Reusable line upper-limit + diagnostics pipeline.

Configurable per target/proposal. Invoke via:
  python -u line_pipeline.py --proposal 2022.1.01160.S --target G189.0307+00.7821 \
    --vlsr 11.5 --distance-kpc 2.0

Outputs land under:
  analysis_products/{target}/{proposal}/
    continuum_sources.csv
    line_measurements.csv
    overview_continuum.png
    overview_lines_topN.png
    source_NN/
      *.spec.npz
      {line}_mom0.fits / _mom1.fits
      {line}_diagnostic.png       # 4-panel per line
      spectrum_full.png            # labelled wide spectrum
      salt_stack.png               # if >=2 NaCl or KCl in band
"""
import argparse
import re
import json
import warnings
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats, mad_std
from astropy.table import Table
from astropy.visualization import simple_norm
from spectral_cube import SpectralCube
from scipy import ndimage as ndi
from photutils.detection import find_peaks
from matplotlib.patches import Circle, Ellipse

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")

VWIN_KMS = 30.0
LINE_MASK_KMS = 25.0
EXT_NBEAM = 3.0
PEAK_SIGMA = 20.0
C_KMS = 299792.458

# ----- Line lists -----
def load_salt_lines(species, ipac_path, fmin=200.0, fmax=280.0, vu_max=2):
    t = Table.read(ipac_path, format="ascii.ipac")
    out = []
    for r in t:
        f = float(r["Frequency"])
        if not (fmin <= f <= fmax): continue
        vu, vl, ju, jl = int(r["vu"]), int(r["vl"]), int(r["Ju"]), int(r["Jl"])
        if vu != vl or (ju - jl) != 1: continue
        if vu > vu_max: continue
        out.append({"name": f"{species}_v{vu}_J{ju}-{jl}",
                    "rest_GHz": f, "Eu_K": float(r["E_U"]),
                    "group": species, "vu": vu, "ju": ju, "jl": jl})
    return out


def build_line_list(fmin=80.0, fmax=500.0):
    """All candidate lines spanning ALMA bands 3-8. Filtered to in-band at runtime."""
    lines = []
    lines += load_salt_lines("NaCl", "/orange/adamginsburg/salt/salt_data/23Na-35Cl_rotational_transitions.ipac", fmin, fmax)
    lines += load_salt_lines("KCl",  "/orange/adamginsburg/salt/salt_data/39K-35Cl_rotational_transitions.ipac",  fmin, fmax)
    # H2O
    lines += [
        # 232.68670 = H2O v2=1 5(5,0)-6(4,3), Eu=3461.9 K (vibrationally
        # excited = genuinely hot water). Name kept as a legacy identifier.
        {"name":"H2O_5_15-4_22_232",      "rest_GHz":232.68670,"Eu_K":3462,"group":"H2O"},
        # H2O_v2_3_13-2_20_232 (232.9366) dropped: always blended with the
        # neighboring CH3OH line, so it never gives a clean H2O measurement.
        {"name":"H2O_5_24-4_31_321",      "rest_GHz":321.22568,"Eu_K":1862,"group":"H2O"},
    ]
    # SO / SO2 / SiO / SiS in 200-360 GHz commonly bright
    other = [
        ("SO_5_6-4_5",     219.94944, 35),
        ("SO_6_5-5_4",     251.82577, 50),
        ("SO_6_7-5_6",     261.84370, 47),
        ("SO_8_7-7_6",     340.71416, 79),
        ("SO_9_8-8_7",     346.52848, 78),
        ("SO2_5_2_4-4_1_3",         241.61580, 24),
        ("SO2_14_4_10-15_3_13",     254.28080, 159),
        ("SO2_18_3_15-18_2_16",     240.94264, 198),
        ("SO2_4_2_2-3_1_3",         235.15170, 19),
        ("SO2_5_3_3-4_2_2",         351.87320, 36),
        ("SiO_v0_J=5-4",   217.10498, 31),
        ("SiO_v0_J=6-5",   260.51800, 44),
        ("SiO_v0_J=8-7",   347.33058, 75),
        ("29SiO_v0_J=6-5", 257.25500, 43),
        ("SiO_v1_J=6-5",   258.70744, 1798),
        ("SiS_v0_13-12",   240.71082, 80),
        ("SiS_v0_14-13",   259.21205, 93),
        ("SiS_v0_19-18",   344.78140, 165),
        ("29SiS_v0_14-13", 256.02400, 92),
    ]
    for n, f, eu in other:
        g = re.match(r"(SO2|SO|SiO|SiS|29SiO|29SiS)", n).group(1)
        lines.append({"name":n,"rest_GHz":f,"Eu_K":eu,"group":g})
    # Recomb
    R_INF_C = 3.289841960364e15
    ME_MP = 0.000544617
    ME_MC = 4.5546e-5   # m_e / M(12C): carbon RRLs sit +149.4 km/s from H
    def rrl(n, dn, me_mx=ME_MP):
        return R_INF_C * (1.0/n**2 - 1.0/(n+dn)**2) / (1.0 + me_mx) * 1e-9
    for dn, sym in [(1, "alpha"), (2, "beta"), (3, "gamma"), (4, "delta")]:
        for n in range(20, 80):
            f = rrl(n, dn)
            if fmin <= f <= fmax:
                lines.append({"name": f"H{n}{sym}", "rest_GHz": f, "Eu_K": 0, "group": "RRL"})
    # Carbon RRLs (alpha only): dense-PDR / photoevaporating-surface tracer.
    # Detected without H counterpart toward IRAS 17008-4040 mm1 (C30alpha
    # 232.017 GHz); keep as a separate group so it isn't conflated with the
    # ionized-gas (H) RRL statistics.
    for n in range(25, 40):
        f = rrl(n, 1, ME_MC)
        if fmin <= f <= fmax:
            lines.append({"name": f"C{n}alpha", "rest_GHz": f, "Eu_K": 0,
                          "group": "CRRL"})
    return lines


# ----- Cube utilities -----
_IRAS_RE = re.compile(r"(?:IRAS[ _]*|\bI)([0-9]{5})", re.IGNORECASE)


def _target_handle(target):
    """Compact target handle used as the filename prefix for per-source
    outputs: prefer 'I12345' IRAS shorthand from sources_L4_d2 alt names,
    else a sanitized version of the target name."""
    src_csv = ROOT / "data/sources_L4_d2.csv"
    if src_csv.exists():
        try:
            df = pd.read_csv(src_csv)
            row = df[df["name"] == target]
            if not row.empty and "alma_target_names" in row.columns:
                alt = str(row.iloc[0]["alma_target_names"] or "")
                m = _IRAS_RE.search(alt)
                if m:
                    return f"I{m.group(1)}"
        except (pd.errors.EmptyDataError, KeyError, OSError):
            pass
    # Sanitize bare target name: keep alnum, +, -, .
    return re.sub(r"[^A-Za-z0-9+\-.]", "", target.replace(" ", ""))


def _fits_getheader_retry(path, attempts=6, delay=5.0):
    """fits.getheader with retry on transient Lustre/FUSE I/O errors."""
    import time
    last_err = None
    for i in range(attempts):
        try:
            return fits.getheader(path)
        except (BrokenPipeError, OSError) as e:
            last_err = e
            time.sleep(delay * (i + 1))
    raise last_err


def cube_freq_range(path):
    try:
        h = _fits_getheader_retry(path)
    except (BrokenPipeError, OSError) as e:
        print(f"  ! cube_freq_range({Path(path).name}): {type(e).__name__}; skipping", flush=True)
        return None, None
    crv = h.get("CRVAL3", 0); cd = h.get("CDELT3", 0); n = h.get("NAXIS3", 0); crp = h.get("CRPIX3", 1)
    if n == 0 or cd == 0: return None, None
    fmin = (crv + (1-crp)*cd)*1e-9
    fmax = (crv + (n-crp)*cd)*1e-9
    return min(fmin,fmax), max(fmin,fmax)


def beam_pix(hdr):
    bmaj = hdr.get("BMAJ", 0); bmin = hdr.get("BMIN", 0)
    cd = abs(hdr.get("CDELT1", hdr.get("CD1_1", 1.0)))
    if cd <= 0 or bmaj <= 0: return 1.0
    return (np.pi * bmaj * bmin) / (4 * np.log(2) * cd ** 2)


# ----- Source ID -----
def _sources_from_region(reg_path, cont_data, cont_hdr, cont_sigma, wcs_cont):
    """Parse a DS9 .reg with circle/ellipse entries (icrs) and return a list
    of source dicts compatible with process_source(). Dedupe identical
    centers (≤0.05" apart). Each source gets peak_Jybeam = continuum value
    at its center pixel.

    Regions emitting `ellipse(ra, dec, a", b", pa)` and `circle(ra, dec, r")`
    are accepted. The aperture's semi-major axis is stored as
    ``aper_arcsec`` for later use when ``--use-region-aperture`` is set.
    """
    import re
    ny, nx = cont_data.shape
    out = []
    seen = []
    raw = Path(reg_path).read_text()
    for ln in raw.splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#") or ln in ("icrs", "fk5"):
            continue
        m_e = re.match(r"ellipse\(([\-0-9.]+),\s*([\-0-9.]+),\s*([\d.]+)\"?,\s*([\d.]+)\"?,\s*([\-0-9.]+)\)", ln)
        m_c = re.match(r"circle\(([\-0-9.]+),\s*([\-0-9.]+),\s*([\d.]+)\"?\)", ln)
        if m_e:
            ra, dec, a, b, pa = (float(m_e.group(i)) for i in range(1, 6))
            aper = max(a, b)
            text_match = re.search(r"text=\{([^}]+)\}", ln)
            label = text_match.group(1) if text_match else ""
        elif m_c:
            ra, dec, r = (float(m_c.group(i)) for i in range(1, 4))
            aper = r
            label = ""
        else:
            continue
        # Dedupe centers within 0.05"
        dup = False
        for ra0, dec0 in seen:
            if abs(ra - ra0) * 3600 * np.cos(np.radians(dec)) < 0.05 and abs(dec - dec0) * 3600 < 0.05:
                dup = True; break
        if dup:
            continue
        seen.append((ra, dec))
        # Pixel coords in continuum image
        xp, yp = wcs_cont.world_to_pixel(SkyCoord(ra * u.deg, dec * u.deg))
        xi, yi = int(round(float(xp))), int(round(float(yp)))
        if not (0 <= xi < nx and 0 <= yi < ny):
            print(f"  region ({ra},{dec}) outside continuum image; skip")
            continue
        peak = float(cont_data[yi, xi])
        out.append(dict(id=len(out) + 1, ra_deg=ra, dec_deg=dec,
                        x=xi, y=yi, peak_Jybeam=peak,
                        snr=peak / cont_sigma if cont_sigma > 0 else 0.0,
                        aper_arcsec=aper, label=label))
    return out


def find_continuum_sources(cont_path, peak_sigma=PEAK_SIGMA, min_peaks=1):
    hdu = fits.open(cont_path)[0]
    data = np.squeeze(hdu.data).astype(np.float32)
    hdr = hdu.header
    wcs = WCS(hdr).celestial
    _, _, sigma = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
    print(f"  continuum sigma = {sigma:.4g} {hdr.get('BUNIT','')}", flush=True)
    threshold = peak_sigma * sigma
    peaks = find_peaks(data, threshold=threshold, box_size=11, npeaks=20)
    if peaks is None or len(peaks) == 0:
        print(f"  no >{peak_sigma}-sigma peaks; lowering to 10-sigma")
        peaks = find_peaks(data, threshold=10*sigma, box_size=11, npeaks=20)
    sources = []
    if peaks is not None:
        for i, p in enumerate(peaks):
            x, y = int(p["x_peak"]), int(p["y_peak"])
            sky = wcs.pixel_to_world(x, y)
            sources.append({
                "id": i+1, "x": x, "y": y,
                "ra_deg": float(sky.ra.deg),
                "dec_deg": float(sky.dec.deg),
                "peak_Jybeam": float(data[y, x]),
                "snr": float(data[y, x] / sigma),
            })
    return sources, sigma, hdr, data, wcs


# ----- Cutout + measure -----
def cutout_cube(cube_path, line, ra, dec, vlsr_kms, half_arcsec=2.5):
    try:
        cube = SpectralCube.read(cube_path)
    except Exception:
        return None
    rest = line["rest_GHz"] * u.GHz
    # RRLs are broad; need wider spectral window for noise + line capture
    vwin = 150.0 if line.get("group") in ("RRL", "CRRL") else VWIN_KMS
    try:
        c = cube.with_spectral_unit(u.km/u.s, velocity_convention="radio", rest_value=rest)
        sub = c.spectral_slab((vlsr_kms - vwin)*u.km/u.s, (vlsr_kms + vwin)*u.km/u.s)
        if sub.shape[0] < 4: return None
    except Exception:
        return None
    wcs_c = cube.wcs.celestial
    coord = SkyCoord(ra, dec, unit="deg")
    try:
        xp, yp = wcs_c.world_to_pixel(coord)
        xp, yp = int(round(float(xp))), int(round(float(yp)))
    except Exception:
        return None
    cd = abs(fits.getheader(cube_path).get("CDELT1", 1)) * 3600
    half = max(8, int(round(half_arcsec / cd)))
    _, ny, nx = sub.shape
    y0, y1 = max(0, yp-half), min(ny, yp+half)
    x0, x1 = max(0, xp-half), min(nx, xp+half)
    if y1-y0 < 6 or x1-x0 < 6: return None
    return sub[:, y0:y1, x0:x1]


def measure_line(cutout, src_pixel_yx, beam_pix_count, vlsr_kms, group="other",
                 on_kms_default=10.0):
    """Group-aware velocity windows: RRLs get broad on-line (+/-60 km/s)
    + far baseline (>100 km/s); other species get +/-10 / >20."""
    data = np.asarray(cutout.unmasked_data[:, :, :].value, dtype=np.float32)
    v = cutout.spectral_axis.to(u.km/u.s).value
    py, px = src_pixel_yx
    py, px = max(0, min(data.shape[1]-1, py)), max(0, min(data.shape[2]-1, px))
    r_pix = max(2.0, np.sqrt(beam_pix_count / np.pi))
    yy, xx = np.indices(data.shape[1:])
    aper = (yy - py)**2 + (xx - px)**2 <= r_pix**2
    spec = np.nanmean(data[:, aper], axis=1)
    # Group-aware windows (off_kms scales with on_kms so the off-line baseline
    # is always 2x the on-line half-window)
    if group in ("RRL", "CRRL"):
        on_kms, off_kms = 60.0, 100.0
    else:
        on_kms = on_kms_default
        off_kms = max(2.0 * on_kms, 20.0)
    off = np.abs(v - vlsr_kms) > off_kms
    if off.sum() < 4: off = np.abs(v - vlsr_kms) > off_kms * 0.6
    if off.sum() < 4: off = np.abs(v - vlsr_kms) > 10.0
    cont = np.nanmedian(spec[off])
    spec_cs = spec - cont
    sigma = mad_std(spec_cs[off]) if off.sum() > 4 else mad_std(spec_cs)
    if not np.isfinite(sigma) or sigma == 0: sigma = np.nanstd(spec_cs)
    on = np.abs(v - vlsr_kms) <= on_kms
    if on.sum() < 2: on = np.abs(v - vlsr_kms) <= VWIN_KMS
    on_data = spec_cs[on] if on.any() else np.array([])
    on_finite = on.any() and np.isfinite(on_data).any()
    peak = float(np.nanmax(on_data)) if on_finite else np.nan
    peak_v = float(v[on][np.nanargmax(on_data)]) if on_finite else np.nan
    dv = float(np.abs(np.median(np.diff(v))))
    integ = float(np.nansum(spec_cs[on]) * dv)
    snr = peak / sigma if sigma > 0 else np.nan
    mom0 = np.nansum((data - cont) * on[:, None, None], axis=0) * dv
    sig_mom0 = mad_std(mom0)
    extended = bool(((mom0 / max(sig_mom0,1e-12)) > 5).sum() / max(beam_pix_count,1) >= EXT_NBEAM)
    return dict(spec=spec_cs, vaxis=v, peak=peak, peak_v=peak_v, sigma=sigma,
                snr=snr, integ=integ, mom0=mom0, sigma_mom0=sig_mom0,
                continuum=cont, extended=extended, raw_data=data,
                aper=aper, beam_pix=beam_pix_count, on_kms=on_kms)


def make_mom1(cutout, vlsr_kms, mom0_thresh_sigma=5.0, on_kms=15.0):
    data = np.asarray(cutout.unmasked_data[:, :, :].value, dtype=np.float32)
    v = cutout.spectral_axis.to(u.km/u.s).value
    on = np.abs(v - vlsr_kms) <= on_kms
    sub = data[on]; vsub = v[on]
    with np.errstate(invalid='ignore', divide='ignore'):
        num = np.nansum(sub * vsub[:, None, None], axis=0)
        den = np.nansum(sub, axis=0)
        m1 = num / den
        m0 = np.nansum(sub, axis=0) * float(np.abs(np.median(np.diff(vsub))))
    sig_m0 = mad_std(m0)
    m1[m0 < mom0_thresh_sigma * sig_m0] = np.nan
    return m1, m0, sig_m0


# ----- Diagnostic plots -----
def plot_continuum_overview(cont_data, cont_hdr, sigma, sources, out_path):
    wcs = WCS(cont_hdr).celestial
    fig = plt.figure(figsize=(10, 9))
    ax = plt.subplot(projection=wcs)
    norm = simple_norm(cont_data, stretch="asinh", min_cut=-3*sigma, max_cut=np.nanpercentile(cont_data, 99.5))
    im = ax.imshow(cont_data, cmap="gray_r", origin="lower", norm=norm)
    # mark sources
    for s in sources:
        circle = Circle((s["x"], s["y"]), radius=8, fill=False, edgecolor="red", lw=1)
        ax.add_patch(circle)
        ax.text(s["x"]+10, s["y"]+10, str(s["id"]), color="red", fontsize=8)
    ax.set_xlabel("RA")
    ax.set_ylabel("Dec")
    ax.set_title(f"Continuum + sources ({len(sources)} >= {PEAK_SIGMA}σ; σ={sigma:.3g})")
    plt.colorbar(im, ax=ax, label=cont_hdr.get("BUNIT", "Jy/beam"))
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def _max_grad_pa_pix(mom1_arr, sy, sx, half=8):
    """Estimate PA (deg, E of N) of max mom1 gradient at (sy, sx)."""
    y0 = max(0, sy - half); y1 = min(mom1_arr.shape[0], sy + half + 1)
    x0 = max(0, sx - half); x1 = min(mom1_arr.shape[1], sx + half + 1)
    sub = mom1_arr[y0:y1, x0:x1]
    if not np.isfinite(sub).any():
        return 0.0
    gy, gx = np.gradient(sub)
    cy, cx = sy - y0, sx - x0
    cy = min(cy, gy.shape[0] - 1); cx = min(cx, gx.shape[1] - 1)
    pa = np.degrees(np.arctan2(gx[cy, cx], gy[cy, cx]))
    return float(pa) if np.isfinite(pa) else 0.0


def plot_line_diagnostic(measure, line, sid, sdir, vlsr_kms, cutout):
    """6-panel: peak | mom0 | mom1+slit | spectrum (w/ contaminants) | PV | info."""
    m = measure
    cube_arr = m["raw_data"]
    v = m["vaxis"]
    mom0 = m["mom0"]
    spec = m["spec"]
    sigma = m["sigma"]
    on_kms = m.get("on_kms", 10.0)
    # compute masked mom1
    sig_m0 = mad_std(mom0) if np.isfinite(mom0).any() else 1.0
    on_mask_v = np.abs(v - vlsr_kms) <= on_kms
    if on_mask_v.any():
        sub = cube_arr[on_mask_v]
        vsub = v[on_mask_v]
        with np.errstate(invalid='ignore', divide='ignore'):
            num = np.nansum(sub * vsub[:, None, None], axis=0)
            den = np.nansum(sub, axis=0)
            mom1 = num / den
        mom1[mom0 < 5*sig_m0] = np.nan
    else:
        mom1 = np.full_like(mom0, np.nan, dtype=np.float32)
    fig, axs = plt.subplots(1, 6, figsize=(24, 4))
    # Compute slit PA from mom1 at source center (cutout-center).
    ny0, nx0 = mom0.shape
    src_yc, src_xc = ny0 // 2, nx0 // 2
    # Crop to the detected signal +-5 px, but never wider than +-3
    # synthesized beams around the SOURCE center (extended contaminants
    # were blowing up signal-based crops and hiding the compact target).
    beam_fwhm_pix = float(np.sqrt(m.get("beam_pix", 25.0) * 4 * np.log(2) / np.pi))
    cap = max(8, int(round(3.0 * beam_fwhm_pix)))
    sig_m0_b = mad_std(mom0, ignore_nan=True) if np.isfinite(mom0).any() else 1.0
    sig_mask = (mom0 > 3.0 * sig_m0_b) & np.isfinite(mom0)
    # only signal within the 3-beam cap participates in the tight crop
    capmask = np.zeros_like(sig_mask)
    capmask[max(0, src_yc - cap):src_yc + cap + 1,
            max(0, src_xc - cap):src_xc + cap + 1] = True
    sig_mask &= capmask
    if sig_mask.sum() >= 3:
        ys_s, xs_s = np.where(sig_mask)
        ymin = max(0, ys_s.min() - 5); ymax = min(ny0, ys_s.max() + 6)
        xmin = max(0, xs_s.min() - 5); xmax = min(nx0, xs_s.max() + 6)
    else:
        half = cap
        ymin = max(0, src_yc - half); ymax = min(ny0, src_yc + half + 1)
        xmin = max(0, src_xc - half); xmax = min(nx0, src_xc + half + 1)
    # arcsec/pix for axis labels
    try:
        cdelt_arcsec_x = abs(cutout.wcs.celestial.wcs.cdelt[0]) * 3600.0
        cdelt_arcsec_y = abs(cutout.wcs.celestial.wcs.cdelt[1]) * 3600.0
    except (AttributeError, IndexError):
        cdelt_arcsec_x = cdelt_arcsec_y = 1.0
    def _arcsec_extent():
        # Center on the source pixel; return (left, right, bottom, top)
        return ((xmin - src_xc) * cdelt_arcsec_x,
                (xmax - src_xc) * cdelt_arcsec_x,
                (ymin - src_yc) * cdelt_arcsec_y,
                (ymax - src_yc) * cdelt_arcsec_y)
    extent = _arcsec_extent()
    # panel 1: peak-channel image
    on = np.abs(v - vlsr_kms) <= 10.0
    if on.any():
        per_chan_max = np.nanmax(cube_arr[on], axis=(1, 2))
        if np.isfinite(per_chan_max).any():
            peak_ch = cube_arr[on][np.nanargmax(per_chan_max)]
        else:
            peak_ch = cube_arr.mean(axis=0)
    else:
        peak_ch = cube_arr.mean(axis=0)
    vmax = np.nanpercentile(peak_ch, 99) if np.isfinite(peak_ch).any() else 1
    vmin = np.nanpercentile(peak_ch, 1)
    peak_crop = peak_ch[ymin:ymax, xmin:xmax]
    axs[0].imshow(peak_crop, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax,
                   extent=extent)
    axs[0].set_title("peak channel")
    axs[0].set_xlabel("$\\Delta$RA (\")"); axs[0].set_ylabel("$\\Delta$Dec (\")")
    axs[0].tick_params(labelsize=8)
    # panel 2: mom0 (tight crop)
    mvmax = max(abs(np.nanpercentile(mom0, 2)), abs(np.nanpercentile(mom0, 98)))
    mom0_crop = mom0[ymin:ymax, xmin:xmax]
    # PuOr (not red/blue): red-blue is reserved for velocity (mom1 panel)
    axs[1].imshow(mom0_crop, origin="lower", cmap="PuOr_r", vmin=-mvmax, vmax=mvmax,
                   extent=extent)
    axs[1].set_title(f"mom0 (±{on_kms:.0f} km/s)\nσ={m['sigma_mom0']:.2g}")
    axs[1].set_xlabel("$\\Delta$RA (\")")
    axs[1].tick_params(labelsize=8)
    # panel 3: masked mom1 + dashed PV slit overlay (tight crop, arcsec axes)
    slit_pa = _max_grad_pa_pix(mom1, src_yc, src_xc) if np.isfinite(mom1).any() else 0.0
    if np.isfinite(mom1).any():
        mom1_crop = mom1[ymin:ymax, xmin:xmax]
        v1, v2 = np.nanpercentile(mom1_crop[np.isfinite(mom1_crop)], [5, 95]) \
                  if np.isfinite(mom1_crop).any() else (-10, 10)
        axs[2].imshow(mom1_crop, origin="lower", cmap="RdBu_r", vmin=v1, vmax=v2,
                       extent=extent)
        axs[2].set_title(f"mom1 (mask mom0>5σ) | PV PA={slit_pa:+.0f}°")
        # Dashed slit line in arcsec coords through source center (0,0).
        half_len_arcsec = max(abs(extent[0]), abs(extent[1]),
                                abs(extent[2]), abs(extent[3])) * 0.9
        pa_rad = np.radians(slit_pa)
        dx_a = np.sin(pa_rad); dy_a = np.cos(pa_rad)
        axs[2].plot([-half_len_arcsec * dx_a, half_len_arcsec * dx_a],
                     [-half_len_arcsec * dy_a, half_len_arcsec * dy_a],
                     "k--", lw=0.8, alpha=0.85)
    else:
        axs[2].text(0.5, 0.5, "no mom1\n(no 5σ pixels)", ha="center", va="center",
                     transform=axs[2].transAxes)
        axs[2].set_title("mom1")
    axs[2].set_xlabel("$\\Delta$RA (\")")
    axs[2].tick_params(labelsize=8)
    # panel 4: spectrum + contaminant labels at observed freq → v conversion
    axs[3].plot(v, spec, "k-", drawstyle="steps-mid", lw=1)
    axs[3].axhline(0, color="grey", lw=0.5)
    axs[3].axhline(sigma, color="red", lw=0.5, ls=":")
    axs[3].axhline(-sigma, color="red", lw=0.5, ls=":")
    axs[3].axhline(3*sigma, color="red", lw=0.5, ls="--")
    axs[3].axvline(vlsr_kms, color="blue", lw=0.5, ls="--")
    # shaded integration bounds
    axs[3].axvspan(vlsr_kms - on_kms, vlsr_kms + on_kms, color="yellow", alpha=0.15)
    axs[3].set_xlabel("v (km/s)")
    axs[3].set_ylabel("intensity")
    axs[3].set_title(f"spec @1-beam\nSNR={m['snr']:.1f}")
    # Overlay contaminant labels: convert each lineid_style rest-line to a
    # velocity-axis position based on THIS line's rest frequency, so that
    # neighboring transitions (NaCl/H2O/ISM/CH3OCHO/CH3OH/RRL) appear at their
    # apparent Doppler velocity offset from the line being plotted.
    try:
        from analysis.lineid_style import (NACL_KCL, WATER, ISM, RRLS,
                                            CH3OCHO_FREQS, CH3OH_FREQS,
                                            COLOR_MAP)
        _, ytop = axs[3].get_ylim()
        ybot, _ = axs[3].get_ylim()
        rest_GHz = float(line["rest_GHz"])
        v_lo, v_hi = float(np.nanmin(v)), float(np.nanmax(v))
        # Show every catalog line whose Doppler-shifted velocity (in this
        # line's frame, accounting for source v_LSR) falls inside the v window
        # of the panel. v_other = (1 - f_other/f_this) * c + v_LSR.
        for (f_other, lab, cls) in (list(NACL_KCL) + list(WATER) + list(ISM)):
            v_other = (1.0 - f_other / rest_GHz) * C_KMS + vlsr_kms
            if not (v_lo <= v_other <= v_hi):
                continue
            col = COLOR_MAP.get(cls, "gray")
            axs[3].axvline(v_other, color=col, ls=":", lw=0.5, alpha=0.7)
            axs[3].text(v_other, ytop, lab, rotation=90, fontsize=6, color=col,
                         ha="right", va="top", alpha=0.85)
        for lab, f_other in RRLS.items():
            v_other = (1.0 - f_other / rest_GHz) * C_KMS + vlsr_kms
            if not (v_lo <= v_other <= v_hi):
                continue
            col = COLOR_MAP.get("rrl", "magenta")
            axs[3].axvline(v_other, color=col, ls=":", lw=0.5, alpha=0.7)
            axs[3].text(v_other, ytop, lab, rotation=90, fontsize=6, color=col,
                         ha="right", va="top", alpha=0.85)
        # CH3OCHO + CH3OH arrows
        yspan = ytop - ybot
        for freqs, arr_cls in ((CH3OCHO_FREQS, "ch3ocho"),
                                (CH3OH_FREQS, "ch3oh")):
            arr_col = COLOR_MAP[arr_cls]
            for f_other in freqs:
                v_other = (1.0 - f_other / rest_GHz) * C_KMS + vlsr_kms
                if not (v_lo <= v_other <= v_hi):
                    continue
                axs[3].annotate("",
                                 xy=(v_other, ytop - 0.18 * yspan),
                                 xytext=(v_other, ytop - 0.02 * yspan),
                                 arrowprops=dict(arrowstyle="-|>",
                                                  color=arr_col, lw=0.6, alpha=0.7))
    except (ImportError, AttributeError, KeyError):
        pass
    # panel 5: PV diagram across max-gradient slit, this line's velocity window
    try:
        from pvextractor import PathFromCenter, extract_pv_slice
        from astropy.coordinates import SkyCoord
        # Use the cutout's center coord as the slit center
        wcs_c = cutout.wcs.celestial
        ra_c, dec_c = wcs_c.pixel_to_world_values(src_xc, src_yc)
        scoord = SkyCoord(ra_c * u.deg, dec_c * u.deg)
        # length: ~4" in arcsec
        cdelt_arcsec = abs(cutout.wcs.celestial.wcs.cdelt[0]) * 3600.0
        slit_len = max(2.0, min(6.0, ny0 * cdelt_arcsec * 0.8))
        path = PathFromCenter(center=scoord,
                              length=slit_len * u.arcsec,
                              angle=slit_pa * u.deg, sample=21)
        # Restrict to ±3×on_kms velocity window for clarity
        cube_v = cutout.with_spectral_unit(u.km / u.s, velocity_convention="radio")
        vmin = (vlsr_kms - 3 * on_kms) * u.km / u.s
        vmax = (vlsr_kms + 3 * on_kms) * u.km / u.s
        sub_v = cube_v.spectral_slab(vmin, vmax)
        pv = extract_pv_slice(sub_v, path)
        pv_data = pv.data
        # Crop offset to ~40 px window so square pixels stay visible
        ny_pv, nx_pv = pv_data.shape
        cx_pv = nx_pv // 2
        crop_h = min(20, cx_pv)
        pv_c = pv_data[:, max(0, cx_pv - crop_h):min(nx_pv, cx_pv + crop_h)]
        pv_sig = float(mad_std(pv_c, ignore_nan=True))
        pv_med = float(np.nanmedian(pv_c))
        # shared PV style: white->black over -3..+3 sigma, black->orange above
        import matplotlib.colors as _mcolors
        _pvmin = pv_med - 3 * pv_sig
        _pvmax = max(float(np.nanmax(pv_c)), pv_med + 6 * pv_sig)
        _fb = min(0.95, max(0.05, (pv_med + 3 * pv_sig - _pvmin) / (_pvmax - _pvmin)))
        _cmap_hc = _mcolors.LinearSegmentedColormap.from_list(
            "pv_hc", [(0.0, "white"), (_fb, "black"), (1.0, "orange")])
        # physical axes: spatial offset (arcsec) x velocity (km/s), plus a
        # secondary frequency axis (channels/pixels are anti-useful)
        c_kms = 299792.458
        lo_i = max(0, cx_pv - crop_h)
        hi_i = min(nx_pv, cx_pv + crop_h)
        cun = str(pv.header.get("CUNIT1", "deg")).lower()
        off_scale = abs(float(pv.header.get("CDELT1", cdelt_arcsec / 3600.0)))
        if cun.startswith("deg"):
            off_scale *= 3600.0
        elif cun.startswith("arcmin"):
            off_scale *= 60.0
        x0 = (lo_i - cx_pv) * off_scale
        x1 = (hi_i - 1 - cx_pv) * off_scale
        vrel = sub_v.spectral_axis.to(u.km / u.s).value - vlsr_kms
        y0, y1 = float(vrel[0]), float(vrel[-1])
        axs[4].imshow(pv_c, origin="lower", aspect="auto",
                       cmap=_cmap_hc, vmin=_pvmin, vmax=_pvmax,
                       extent=[x0, x1, y0, y1])
        axs[4].axvline(0.0, color="0.4", linestyle="--",
                        linewidth=0.8, alpha=0.8)
        axs[4].set_title(f"PV (line slit)\nlen={slit_len:.1f}\"")
        axs[4].set_xlabel("offset (\")")
        axs[4].set_ylabel("v - v$_{lsr}$ (km/s)")
        _rest = float(line["rest_GHz"])
        _v2f = lambda vr: _rest * (1.0 - (vr + vlsr_kms) / c_kms)
        _f2v = lambda f: (1.0 - f / _rest) * c_kms - vlsr_kms
        _sec = axs[4].secondary_yaxis("right", functions=(_v2f, _f2v))
        _sec.set_ylabel("sky freq (GHz)")
    except (ImportError, ValueError, RuntimeError, OSError) as _pverr:
        # Keep the PV panel visible (framed) even on failure so every
        # diagnostic has a PV slot, per the requirement that all diagnostics
        # (detection or not) include a PV panel.
        axs[4].text(0.5, 0.5, f"PV unavailable:\n{_pverr}",
                     ha="center", va="center",
                     transform=axs[4].transAxes, fontsize=8)
        axs[4].set_title("PV (line slit)")
        axs[4].set_xlabel("offset (\")")
        axs[4].set_ylabel("v - v$_{lsr}$ (km/s)")
    # panel 6: text summary
    axs[5].axis("off")
    info = (f"line: {line['name']}\n"
            f"rest: {line['rest_GHz']:.5f} GHz\n"
            f"Eu: {line['Eu_K']:.1f} K\n"
            f"group: {line['group']}\n"
            f"peak: {m['peak']:.3g}\n"
            f"σ (1-beam): {m['sigma']:.3g}\n"
            f"SNR: {m['snr']:.2f}\n"
            f"peak_v: {m['peak_v']:.2f} km/s\n"
            f"integ: {m['integ']:.3g}\n"
            f"extended: {m['extended']}\n"
            f"vlsr (ref): {vlsr_kms:.1f}\n")
    axs[5].text(0.02, 0.98, info, transform=axs[5].transAxes, va="top",
                family="monospace", fontsize=9)
    # Reference position in suptitle (cutout center sky coord)
    try:
        ref_ra, ref_dec = cutout.wcs.celestial.pixel_to_world_values(src_xc, src_yc)
        coord_str = f"(α={float(ref_ra):.5f}°, δ={float(ref_dec):+.5f}°)"
    except (AttributeError, ValueError):
        coord_str = ""
    fig.suptitle(f"src {sid} | {line['name']}  {coord_str}", fontsize=11)
    fig.tight_layout()
    fig.savefig(sdir / f"{line['name']}_diagnostic.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def plot_full_spectrum(sid, sdir, rows, lines_by_name, vlsr_kms):
    """One labelled spectrum stitching all in-band lines (per-line panels)."""
    inband = [r for r in rows if r.get("in_band") and "peak_Kkms_or_unit" in r]
    if not inband:
        return
    n = len(inband)
    ncol = 3
    nrow = (n + ncol - 1) // ncol
    fig, axs = plt.subplots(nrow, ncol, figsize=(4*ncol, 2.5*nrow))
    axs = np.atleast_2d(axs).reshape(nrow, ncol)
    for i, r in enumerate(inband):
        ax = axs[i//ncol, i%ncol]
        sf = sdir / f"{r['line']}.spec.npz"
        if not sf.exists():
            ax.set_visible(False); continue
        arr = np.load(sf)
        v = arr["vaxis"]; spec = arr["spec"]; sg = float(arr["sigma"])
        ax.plot(v, spec, "k-", drawstyle="steps-mid", lw=0.8)
        ax.axhline(0, color="grey", lw=0.4)
        ax.axhline(3*sg, color="red", lw=0.4, ls=":")
        ax.axhline(-3*sg, color="red", lw=0.4, ls=":")
        ax.axvline(vlsr_kms, color="blue", lw=0.4, ls="--")
        line = lines_by_name[r["line"]]
        snr = r.get("snr", np.nan)
        snr_s = f"{snr:.1f}" if np.isfinite(snr) else "n/a"
        ax.set_title(f"{r['line']} ({line['group']}) {snr_s}σ", fontsize=8)
        ax.tick_params(labelsize=7)
    # blank unused
    for j in range(len(inband), nrow*ncol):
        axs[j//ncol, j%ncol].set_visible(False)
    fig.suptitle(f"Source {sid} - per-line aperture-mean spectra (NOT stacked; vlsr={vlsr_kms:.1f} km/s)",
                  fontsize=11)
    fig.tight_layout()
    fig.savefig(sdir / "spectrum_full.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def plot_per_group_stack(sid, sdir, rows, lines_by_name, vlsr_kms):
    """Velocity-shifted-and-stacked spectra per species group (RRL, SO2, NaCl, etc.).
    Single figure with one panel per group."""
    inband = [r for r in rows if r.get("in_band") and "peak_Kkms_or_unit" in r]
    if not inband:
        return
    by_group = {}
    for r in inband:
        g = r["group"]
        sf = sdir / f"{r['line']}.spec.npz"
        if not sf.exists():
            continue
        arr = np.load(sf)
        v = arr["vaxis"] - vlsr_kms
        order = np.argsort(v)
        v_sorted = v[order]
        s_sorted = arr["spec"][order]
        sigma = float(arr["sigma"])
        by_group.setdefault(g, []).append((v_sorted, s_sorted, sigma))
    v_grid = np.arange(-VWIN_KMS, VWIN_KMS + 0.5, 1.0)
    stacks = {}
    for g, entries in by_group.items():
        interps, sigs = [], []
        for v_s, s_s, sig in entries:
            interp = np.interp(v_grid, v_s, s_s, left=np.nan, right=np.nan)
            if np.isfinite(interp).any():
                interps.append(interp); sigs.append(sig)
        if not interps:
            continue
        w = 1.0 / np.array(sigs) ** 2
        stack = np.nansum(np.array(interps) * w[:, None], axis=0) / np.nansum(w)
        sig_stack = 1.0 / np.sqrt(np.nansum(w))
        # Skip degenerate
        if not (np.isfinite(stack) & (stack != 0)).any() or np.nanstd(stack) == 0:
            continue
        stacks[g] = (stack, sig_stack, len(interps))
    if not stacks:
        return
    n = len(stacks)
    ncol = 2
    nrow = (n + ncol - 1) // ncol
    fig, axs = plt.subplots(nrow, ncol, figsize=(5.5*ncol, 2.5*nrow), sharex=True)
    axs = np.atleast_2d(axs).reshape(nrow, ncol)
    for i, (g, (stack, sig, n_lines)) in enumerate(sorted(stacks.items())):
        ax = axs[i//ncol, i%ncol]
        ax.plot(v_grid, stack, "k-", drawstyle="steps-mid", lw=1)
        ax.axhline(0, color="grey", lw=0.4)
        ax.axhline(sig, color="red", lw=0.4, ls=":")
        ax.axhline(-sig, color="red", lw=0.4, ls=":")
        ax.axhline(3*sig, color="red", lw=0.4, ls="--")
        ax.axvline(0, color="blue", lw=0.4, ls="--")
        ax.set_title(f"{g} stack (N={n_lines}, σ={sig:.2g})", fontsize=9)
        ax.set_ylabel("intensity")
    for j in range(len(stacks), nrow*ncol):
        axs[j//ncol, j%ncol].set_visible(False)
    for ax in axs[-1, :]:
        ax.set_xlabel("v - v_lsr (km/s)")
    fig.suptitle(f"Source {sid} - per-group velocity-shifted-and-stacked spectra (vlsr={vlsr_kms:.1f})",
                  fontsize=11)
    fig.tight_layout()
    fig.savefig(sdir / "spectrum_per_group_stack.png", dpi=120, bbox_inches="tight")
    plt.close(fig)


def _plot_one_stack(ax, d, title):
    """Single-panel stack drawing with included-line label list."""
    ax.plot(d["vaxis"], d["spec"], "k-", drawstyle="steps-mid", lw=1)
    ax.axhline(0, color="grey", lw=0.5)
    ax.axhline(d["sigma"], color="red", lw=0.5, ls=":",
                label=f"1σ={d['sigma']:.2g}")
    ax.axhline(-d["sigma"], color="red", lw=0.5, ls=":")
    ax.axhline(3 * d["sigma"], color="red", lw=0.5, ls="--",
                 label=f"3σ={3*d['sigma']:.2g}")
    ax.axvline(0, color="blue", lw=0.5, ls="--", label="v_lsr")
    n = len(d.get("lines_used", []))
    ax.set_ylabel(f"{title} stack\n(N={n}, T_K)")
    # Annotate the line names actually included in the stack
    if "lines_used" in d and d["lines_used"]:
        # Wrap the line list at every 5 entries
        lines_used = list(d["lines_used"])
        chunks = [", ".join(lines_used[i:i+5]) for i in range(0, len(lines_used), 5)]
        ax.text(0.01, 0.98, "lines stacked:\n" + "\n".join(chunks),
                 transform=ax.transAxes, ha="left", va="top", fontsize=7,
                 color="navy",
                 bbox=dict(facecolor="white", alpha=0.7, edgecolor="none",
                            pad=2))
    ax.legend(loc="upper right", fontsize=7)


def plot_salt_stack(stacks, sdir, sid):
    """Write three stack PNGs (when populated):

      nacl_stack.png      : single-panel NaCl-v0/v1/v2-35Cl
      kcl_stack.png       : single-panel KCl-v0/v1/v2-35Cl
      naclkcl_stack.png   : two-panel joint NaCl + KCl

    Each panel labels the line names used in the stack.
    """
    if not stacks:
        return
    # NaCl-only
    if "NaCl" in stacks:
        fig, ax = plt.subplots(figsize=(10, 3))
        _plot_one_stack(ax, stacks["NaCl"], "NaCl 35Cl v=0/1/2")
        ax.set_xlabel("v - v_lsr (km/s)")
        fig.suptitle(f"Source {sid}: NaCl stacked spectrum")
        fig.tight_layout()
        fig.savefig(sdir / "nacl_stack.png", dpi=130, bbox_inches="tight")
        plt.close(fig)
    # KCl-only
    if "KCl" in stacks:
        fig, ax = plt.subplots(figsize=(10, 3))
        _plot_one_stack(ax, stacks["KCl"], "KCl 35Cl v=0/1/2")
        ax.set_xlabel("v - v_lsr (km/s)")
        fig.suptitle(f"Source {sid}: KCl stacked spectrum")
        fig.tight_layout()
        fig.savefig(sdir / "kcl_stack.png", dpi=130, bbox_inches="tight")
        plt.close(fig)
    # Joint NaCl + KCl: two-panel
    species_present = [s for s in ("NaCl", "KCl") if s in stacks]
    if len(species_present) >= 1:
        fig, axs = plt.subplots(len(species_present), 1,
                                  figsize=(10, 3 * len(species_present)),
                                  sharex=True)
        if len(species_present) == 1:
            axs = [axs]
        for ax, sp in zip(axs, species_present):
            _plot_one_stack(ax, stacks[sp], f"{sp} 35Cl v=0/1/2")
        axs[-1].set_xlabel("v - v_lsr (km/s)")
        fig.suptitle(f"Source {sid}: joint NaCl+KCl stacked spectra")
        fig.tight_layout()
        fig.savefig(sdir / "naclkcl_stack.png", dpi=130, bbox_inches="tight")
        plt.close(fig)
    # Also build the JOINT NaCl+KCl single-stack: concat all lines
    d_c = None
    if "NaCl" in stacks and "KCl" in stacks:
        d1, d2 = stacks["NaCl"], stacks["KCl"]
        w1 = 1.0 / d1["sigma"] ** 2
        w2 = 1.0 / d2["sigma"] ** 2
        combined = (d1["spec"] * w1 + d2["spec"] * w2) / (w1 + w2)
        sigma_c = 1.0 / np.sqrt(w1 + w2)
        n_c = len(d1.get("lines_used", [])) + len(d2.get("lines_used", []))
        lines_c = list(d1.get("lines_used", [])) + list(d2.get("lines_used", []))
        d_c = dict(vaxis=d1["vaxis"], spec=combined, sigma=sigma_c,
                    n_lines=n_c, lines_used=lines_c)
        fig, ax = plt.subplots(figsize=(10, 3))
        _plot_one_stack(ax, d_c, "NaCl+KCl 35Cl v=0/1/2 combined")
        ax.set_xlabel("v - v_lsr (km/s)")
        fig.suptitle(f"Source {sid}: joint NaCl+KCl 35Cl-v=0/1/2 stacked spectrum")
        fig.tight_layout()
        fig.savefig(sdir / "naclkcl_combined_stack.png", dpi=130,
                     bbox_inches="tight")
        plt.close(fig)
        np.savez(sdir / "naclkcl_combined_stack.npz",
                  vaxis=d_c["vaxis"], spec=d_c["spec"],
                  sigma=d_c["sigma"], n_lines=d_c["n_lines"])

    # 3-panel salt_stack.png: joint(NaCl+KCl) on top, then NaCl, then KCl.
    # When only one species is present, render only the panels that have data.
    panels_three = []
    if d_c is not None:
        panels_three.append((d_c, "NaCl+KCl combined"))
    if "NaCl" in stacks:
        panels_three.append((stacks["NaCl"], "NaCl 35Cl v=0/1/2"))
    if "KCl" in stacks:
        panels_three.append((stacks["KCl"], "KCl 35Cl v=0/1/2"))
    if panels_three:
        n = len(panels_three)
        fig, axs = plt.subplots(n, 1, figsize=(10, 2.7 * n), sharex=True)
        if n == 1:
            axs = [axs]
        for ax, (d, title) in zip(axs, panels_three):
            _plot_one_stack(ax, d, title)
        axs[-1].set_xlabel("v - v_lsr (km/s)")
        fig.suptitle(f"Source {sid}: salt stacked spectra")
        fig.tight_layout()
        fig.savefig(sdir / "salt_stack.png", dpi=130, bbox_inches="tight")
        plt.close(fig)


# ----- Stacking -----
def stack_salt(rows, sdir, vlsr_kms):
    """Stack the 35Cl isotopologue v=0/1/2 lines (the common, brightest set)
    for NaCl and KCl. Excludes 37Cl/41K/v>=3 lines so the stacks compare
    apples to apples across sources.

    Returns dict {species: {vaxis, spec, sigma, n_lines, lines_used}}."""
    import re as _re
    out = {}
    isotop_re = _re.compile(r"^(NaCl|KCl)_v[012]_")
    for species in ("NaCl", "KCl"):
        sub = [r for r in rows if r.get("group") == species and r.get("in_band")
               and "peak_Kkms_or_unit" in r
               and isotop_re.match(str(r.get("line", "")))]
        if len(sub) < 2:
            continue
        stacks, sigmas, used = [], [], []
        v_grid = np.arange(-VWIN_KMS, VWIN_KMS + 0.5, 1.0)
        for r in sub:
            sf = sdir / f"{r['line']}.spec.npz"
            if not sf.exists():
                continue
            arr = np.load(sf)
            v_shift = arr["vaxis"] - vlsr_kms
            order = np.argsort(v_shift)
            interp = np.interp(v_grid, v_shift[order], arr["spec"][order],
                                left=np.nan, right=np.nan)
            if not np.isfinite(interp).any():
                continue
            stacks.append(interp)
            sigmas.append(float(arr["sigma"]))
            used.append(str(r["line"]))
        if len(stacks) < 2:
            continue
        weights = 1.0 / np.array(sigmas) ** 2
        stack = np.nansum(np.array(stacks) * weights[:, None], axis=0) / np.nansum(weights)
        sigma_stack = 1.0 / np.sqrt(np.nansum(weights))
        finite_nonzero = np.isfinite(stack) & (stack != 0)
        if finite_nonzero.sum() < 5 or np.nanstd(stack) == 0:
            print(f"  ! {species} stack degenerate; not saving")
            continue
        out[species] = dict(vaxis=v_grid, spec=stack, sigma=sigma_stack,
                              n_lines=len(stacks), lines_used=used)
        np.savez(sdir / f"{species}_stack.npz",
                  vaxis=v_grid, spec=stack, sigma=sigma_stack,
                  n_lines=len(stacks), lines_used=np.array(used))
    return out


# ----- Per-source processing -----
def process_source(src, cubes, lines, outroot, vlsr_kms, distance_kpc, on_kms_default=10.0):
    sid = src["id"]
    # Robust filename label: <handle>mm<rank> avoids cross-proposal collisions
    # where per-proposal numeric ids point to different physical sources.
    handle = src.get("handle", f"src{sid:02d}")
    rank = src.get("rank", sid)
    label = f"{handle}mm{rank}"
    sdir = outroot / label
    sdir.mkdir(parents=True, exist_ok=True)
    # Compatibility shim: ensure a legacy source_NN/ symlink still exists so
    # CARTA snippets and figure scripts that hardcode source_NN/ don't break.
    legacy = outroot / f"source_{sid:02d}"
    if not legacy.exists():
        try:
            legacy.symlink_to(label)
        except OSError:
            pass
    print(f"\n=== {label} (id={sid}): ra={src['ra_deg']:.5f} dec={src['dec_deg']:.5f} "
          f"cont={src['peak_Jybeam']:.3g} ({src['snr']:.1f}σ) ===", flush=True)
    rows = []
    lines_by_name = {l["name"]: l for l in lines}
    half_arcsec = max(1.5, 2000.0 / (distance_kpc*1000.0))  # >=1.5" or 2000 AU
    for line in lines:
        rest = line["rest_GHz"]
        f_obs = rest * (1 - vlsr_kms/C_KMS)
        cube_match = None
        for cp in cubes:
            fmin, fmax = cube_freq_range(cp)
            if fmin is None: continue
            if fmin + 0.05 <= f_obs <= fmax - 0.05:
                cube_match = cp; break
        if cube_match is None:
            rows.append({"source": sid, "line": line["name"], "group": line["group"],
                         "rest_GHz": rest, "Eu_K": line["Eu_K"], "in_band": False})
            continue
        cutout = cutout_cube(cube_match, line, src["ra_deg"], src["dec_deg"], vlsr_kms,
                              half_arcsec=half_arcsec)
        if cutout is None:
            rows.append({"source": sid, "line": line["name"], "group": line["group"],
                         "rest_GHz": rest, "Eu_K": line["Eu_K"], "in_band": True,
                         "cutout_fail": True})
            continue
        bp = beam_pix(fits.getheader(cube_match))
        wcs_c = cutout.wcs.celestial
        coord = SkyCoord(src["ra_deg"], src["dec_deg"], unit="deg")
        xp, yp = wcs_c.world_to_pixel(coord)
        peak_yx = (int(round(float(yp))), int(round(float(xp))))
        m = measure_line(cutout, peak_yx, bp, vlsr_kms, group=line["group"],
                          on_kms_default=on_kms_default)
        row = {"source": sid, "line": line["name"], "group": line["group"],
               "rest_GHz": rest, "Eu_K": line["Eu_K"], "in_band": True,
               "peak_Kkms_or_unit": float(m["peak"]),
               "sigma": float(m["sigma"]),
               "snr": float(m["snr"]) if np.isfinite(m["snr"]) else np.nan,
               "integ": float(m["integ"]),
               "peak_v": float(m["peak_v"]) if np.isfinite(m["peak_v"]) else np.nan,
               "extended": bool(m["extended"]),
               "continuum_baseline": float(m["continuum"]),
               "cube": Path(cube_match).name}
        rows.append(row)
        np.savez(sdir / f"{line['name']}.spec.npz", vaxis=m["vaxis"], spec=m["spec"], sigma=m["sigma"])
        # Always make mom0+mom1 FITS + diagnostic PNG
        m1, m0, sig_m0 = make_mom1(cutout, vlsr_kms, mom0_thresh_sigma=3.0,
                                     on_kms=70.0 if line["group"] in ("RRL", "CRRL") else 15.0)
        hdr_m = cutout.wcs.celestial.to_header()
        hdr_m["BUNIT"] = "km/s"
        fits.PrimaryHDU(m1, hdr_m).writeto(sdir / f"{label}_{line['name']}_mom1.fits", overwrite=True)
        hdr0 = hdr_m.copy(); hdr0["BUNIT"] = "Jy/beam km/s"
        fits.PrimaryHDU(m0, hdr0).writeto(sdir / f"{label}_{line['name']}_mom0.fits", overwrite=True)
        try:
            plot_line_diagnostic(m, line, sid, sdir, vlsr_kms, cutout)
        except Exception as e:
            print(f"  ! diagnostic plot fail {line['name']}: {e}")
    # Full source spectrum overview (per-line aperture-mean grid, NOT stacked)
    try:
        plot_full_spectrum(sid, sdir, rows, lines_by_name, vlsr_kms)
    except Exception as e:
        print(f"  ! full spectrum plot fail: {e}")
    # Per-group velocity-shifted-and-stacked spectra
    try:
        plot_per_group_stack(sid, sdir, rows, lines_by_name, vlsr_kms)
    except Exception as e:
        print(f"  ! per-group stack plot fail: {e}")
    return rows, sdir


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proposal", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--vlsr", type=float, required=True)
    ap.add_argument("--distance-kpc", type=float, required=True)
    ap.add_argument("--peak-sigma", type=float, default=PEAK_SIGMA)
    ap.add_argument("--on-kms", type=float, default=10.0,
                    help="non-RRL line integration half-window (km/s); "
                         "RRLs always use 60 km/s")
    ap.add_argument("--sources-region", default=None,
                    help="DS9 .reg with explicit source apertures (icrs); "
                         "circle/ellipse centers become sources, replacing "
                         "the continuum auto-detect. Aperture is still "
                         "the cube beam unless --use-region-aperture is set.")
    ap.add_argument("--use-region-aperture", action="store_true",
                    help="When --sources-region is set, use each region's "
                         "ellipse (semi-major axis) as the aperture radius "
                         "instead of the cube beam.")
    ap.add_argument("--target-tokens", default="",
                    help="comma-separated extra filename substrings that "
                         "identify the right target (e.g. 'I09002-4732' for "
                         "G268.4222-00.8490). Used when the MOUS contains "
                         "multiple targets so the pipeline picks the right "
                         "continuum + cubes.")
    args = ap.parse_args()

    srcdir = ROOT / f"uvdata/{args.proposal}/{args.target}"
    outroot = ROOT / f"analysis_products/{args.target}/{args.proposal}"
    outroot.mkdir(parents=True, exist_ok=True)

    # Filename-target matching to handle multi-source MOUS (e.g., 2022.1.01344.S)
    # ALMA encodes + as p and - as m in filenames.
    target_variants = {args.target,
                       args.target.replace("+", "p").replace("-", "m"),
                       args.target.replace("+", "").replace("-", "")}
    for tok in args.target_tokens.split(","):
        tok = tok.strip()
        if tok:
            target_variants.add(tok)
    def name_matches_target(path):
        nm = path.name.lower()
        return any(v.lower() in nm for v in target_variants)

    all_cubes = sorted(srcdir.glob("member.*sci.*cube*pbcor.fits"))
    cubes_target = [c for c in all_cubes if name_matches_target(c)]
    cubes = cubes_target if cubes_target else all_cubes
    print(f"Cubes: {len(all_cubes)} total, {len(cubes_target)} target-matched -> using {len(cubes)}", flush=True)

    # Preferred: combined-SPW continuum (.cont.*)
    cont_candidates = list(srcdir.glob("member.*sci.spw*_*.cont*tt0.pbcor.fits")) \
                    + list(srcdir.glob("member.*sci.spw*_*.cont*.I.pbcor.fits"))
    # Fallback 1: any single-SPW mfs map
    if not cont_candidates:
        cont_candidates = list(srcdir.glob("member.*sci.spw*.mfs.*.I.pbcor.fits")) \
                        + list(srcdir.glob("member.*sci.spw*.mfs.I.pbcor.fits"))
        if cont_candidates:
            print(f"  no .cont — falling back to mfs map", flush=True)
    cont_target = [c for c in cont_candidates if name_matches_target(c)]
    cont_candidates = cont_target if cont_target else cont_candidates
    # Fallback 2: derive cont from first cube via spectral_cube median
    if not cont_candidates:
        if not cubes:
            raise FileNotFoundError(f"No continuum or cubes in {srcdir}")
        print(f"  no .cont or mfs — deriving cont from cube channels", flush=True)
        from spectral_cube import SpectralCube as _SC
        cube0 = _SC.read(cubes[0])
        # Median of first/last 5% channels
        nch = cube0.shape[0]
        edge = max(2, nch // 20)
        slab1 = cube0[:edge, :, :]
        slab2 = cube0[-edge:, :, :]
        d1 = np.asarray(slab1.unmasked_data[:, :, :].value, dtype=np.float32)
        d2 = np.asarray(slab2.unmasked_data[:, :, :].value, dtype=np.float32)
        cont_arr = np.nanmedian(np.concatenate([d1, d2], axis=0), axis=0)
        cont_hdr0 = fits.getheader(cubes[0])
        derived_path = outroot / f"derived_continuum_from_{Path(cubes[0]).name}"
        # Need 2D header
        from astropy.wcs import WCS as _WCS
        wcs2d = _WCS(cont_hdr0).celestial
        hdr2d = wcs2d.to_header()
        for k in ("BMAJ","BMIN","BPA","BUNIT","BTYPE"):
            if k in cont_hdr0:
                hdr2d[k] = cont_hdr0[k]
        fits.PrimaryHDU(cont_arr, hdr2d).writeto(derived_path, overwrite=True)
        cont_candidates = [derived_path]
    cont_path = cont_candidates[0]
    print(f"Cont: {cont_path.name}")
    print(f"Cubes: {len(cubes)}")

    # Build line list constrained to band
    fmins, fmaxs = [], []
    for cp in cubes:
        fl, fh = cube_freq_range(cp)
        if fl: fmins.append(fl); fmaxs.append(fh)
    fmin_all = min(fmins) - 0.5 if fmins else 80.0
    fmax_all = max(fmaxs) + 0.5 if fmaxs else 500.0
    lines = build_line_list(fmin=fmin_all, fmax=fmax_all)
    print(f"Candidate lines: {len(lines)} (in 80-500 GHz, span {fmin_all:.1f}-{fmax_all:.1f})")

    if args.sources_region:
        sources_auto, cont_sigma, cont_hdr, cont_data, wcs_cont = find_continuum_sources(
            cont_path, args.peak_sigma)
        sources = _sources_from_region(args.sources_region, cont_data, cont_hdr,
                                         cont_sigma, wcs_cont)
        print(f"  {len(sources)} manual-region sources from "
              f"{args.sources_region}", flush=True)
    else:
        sources, cont_sigma, cont_hdr, cont_data, wcs_cont = find_continuum_sources(cont_path, args.peak_sigma)
        print(f"  {len(sources)} sources >= {args.peak_sigma}σ", flush=True)
    # Assign a target-prefixed handle + brightness rank to each source so
    # filenames are unique across proposals. Brightness rank is 1 for the
    # brightest peak_Jybeam, 2 for the next, etc.
    handle = _target_handle(args.target)
    if sources:
        order = sorted(range(len(sources)),
                       key=lambda i: -float(sources[i].get("peak_Jybeam", 0.0)))
        for rank_minus_one, idx in enumerate(order):
            sources[idx]["handle"] = handle
            sources[idx]["rank"] = rank_minus_one + 1
    pd.DataFrame(sources).to_csv(outroot / "continuum_sources.csv", index=False)

    plot_continuum_overview(cont_data, cont_hdr, cont_sigma, sources, outroot / "overview_continuum.png")

    all_rows = []
    # Per-source vlsr overrides (canonical args.vlsr as default; refine per-source from RRL/SO2 peak_v)
    vlsr_per_source = {}
    src_overrides_path = outroot / "vlsr_per_source.csv"
    if src_overrides_path.exists():
        import csv as _csv
        with open(src_overrides_path) as fh:
            for row in _csv.DictReader(fh):
                vlsr_per_source[int(row["source"])] = float(row["vlsr_kms"])
        print(f"Loaded vlsr overrides for {len(vlsr_per_source)} sources from {src_overrides_path.name}")

    for s in sources:
        v_use = vlsr_per_source.get(s["id"], args.vlsr)
        rows, sdir = process_source(s, cubes, lines, outroot, v_use, args.distance_kpc,
                                     on_kms_default=args.on_kms)
        stacks = stack_salt(rows, sdir, args.vlsr)
        plot_salt_stack(stacks, sdir, s["id"])
        all_rows.extend(rows)

    out_csv = outroot / "line_measurements.csv"
    df_all = pd.DataFrame(all_rows)
    if df_all.empty or "in_band" not in df_all.columns:
        df_keep = df_all
    else:
        df_keep = df_all[df_all.in_band == True].drop(columns=["in_band"])
    df_keep.to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv} ({len(df_keep)} rows; dropped {len(df_all)-len(df_keep)} out-of-band)", flush=True)

    # Source region file (DS9 fk5)
    reg_path = outroot / "sources.reg"
    with open(reg_path, "w") as fh:
        fh.write("# Region file format: DS9 version 4.1\n")
        fh.write("global color=red dashlist=8 3 width=1 font=\"helvetica 10 normal\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
        fh.write("fk5\n")
        for s in sources:
            fh.write(f'circle({s["ra_deg"]:.6f},{s["dec_deg"]:.6f},0.5") # text={{src{s["id"]}}}\n')
    print(f"wrote {reg_path}", flush=True)

    # Top-3 mom0 maps by SNR
    if df_keep.empty or "snr" not in df_keep.columns:
        df_in = df_keep
        top3 = df_keep.iloc[0:0]
    else:
        df_in = df_keep.dropna(subset=["snr"])
        top3 = df_in.nlargest(3, "snr")
    top3_paths = []
    for _, r in top3.iterrows():
        sid = int(r['source'])
        sdir = outroot / f"source_{sid:02d}"
        # New naming: source_<NN>_<line>_mom0.fits. Fall back to legacy
        # <line>_mom0.fits for pre-existing analysis runs.
        p_new = sdir / f"source_{sid:02d}_{r['line']}_mom0.fits"
        p_old = sdir / f"{r['line']}_mom0.fits"
        if p_new.exists():
            top3_paths.append(p_new)
        elif p_old.exists():
            top3_paths.append(p_old)
    print(f"top-3 mom0 maps: {[p.name for p in top3_paths]}")

    # CARTA snippet — paths must be relative to CARTA root (/orange/adamginsburg)
    snip_dir = Path.home() / ".carta/config/snippets"
    snip_dir.mkdir(parents=True, exist_ok=True)
    CARTA_PREFIX = "/orange/adamginsburg"

    def _carta(p):
        s = str(p)
        return s[len(CARTA_PREFIX):] if s.startswith(CARTA_PREFIX) else s

    # Consolidated per-target snippet name (no proposal suffix); merge with
    # existing entries so per-proposal runs accumulate into one snippet.
    snip_out = snip_dir / f"{args.target}_2026.json"
    new_lines = [f'await app.openFile("{_carta(cont_path)}")']
    for cube in cubes:
        new_lines.append(f'await app.appendFile("{_carta(cube)}")')
    for p in top3_paths:
        new_lines.append(f'await app.appendFile("{_carta(p)}")')
    reg_dir = str(Path(reg_path).parent) + "/"
    reg_name = Path(reg_path).name
    new_lines.append(f'await app.importRegion("{_carta(reg_dir)}", "{reg_name}", 2)')

    existing_lines = []
    if snip_out.exists():
        try:
            existing = json.loads(snip_out.read_text())
            existing_lines = (existing.get("code", "") or "").split("\n")
        except json.JSONDecodeError:
            existing_lines = []
    # First line must be openFile; everything else can be appendFile. We dedup
    # by file path so reruns don't grow unbounded.
    seen = set()
    out_lines = []
    for line in existing_lines + new_lines:
        line = line.strip()
        if not line:
            continue
        # Convert any leftover openFile-on-non-first to appendFile after we
        # already have at least one open call recorded.
        if out_lines and line.startswith('await app.openFile('):
            line = line.replace('await app.openFile(', 'await app.appendFile(', 1)
        # dedup key: the path inside the call (or full call for importRegion)
        key = line
        if key in seen:
            continue
        seen.add(key)
        out_lines.append(line)
    snippet = {
        "$schema": "https://cartavis.github.io/schemas/snippet_schema_1.json",
        "categories": ["saltsearch2026"],
        "code": "\n".join(out_lines),
        "frontendVersion": "3.0.0",
        "snippetVersion": 1,
    }
    snip_out.write_text(json.dumps(snippet, indent=4))
    print(f"wrote CARTA snippet: {snip_out}")


if __name__ == "__main__":
    main()
