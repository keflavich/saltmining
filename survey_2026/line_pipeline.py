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
        {"name":"H2O_5_15-4_22_232",      "rest_GHz":232.68670,"Eu_K":644,"group":"H2O"},
        {"name":"H2O_v2_3_13-2_20_232",   "rest_GHz":232.93660,"Eu_K":2400,"group":"H2O"},
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
    def rrl(n, dn): return R_INF_C * (1.0/n**2 - 1.0/(n+dn)**2) / (1.0 + ME_MP) * 1e-9
    for dn, sym in [(1, "alpha"), (2, "beta"), (3, "gamma"), (4, "delta")]:
        for n in range(20, 80):
            f = rrl(n, dn)
            if fmin <= f <= fmax:
                lines.append({"name": f"H{n}{sym}", "rest_GHz": f, "Eu_K": 0, "group": "RRL"})
    return lines


# ----- Cube utilities -----
def cube_freq_range(path):
    h = fits.getheader(path)
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
    vwin = 150.0 if line.get("group") == "RRL" else VWIN_KMS
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
    if group == "RRL":
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


def plot_line_diagnostic(measure, line, sid, sdir, vlsr_kms, cutout):
    """5-panel: peak channel | mom0 | masked mom1 | spectrum (w/ integ bounds) | info"""
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
    fig, axs = plt.subplots(1, 5, figsize=(20, 4))
    # 3x crop on the image panels (panels 0, 1, 2)
    ZOOM = 3.0
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
    axs[0].imshow(peak_ch, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax)
    axs[0].set_title("peak channel")
    axs[0].set_xticks([]); axs[0].set_yticks([])
    # panel 2: mom0
    mvmax = max(abs(np.nanpercentile(mom0, 2)), abs(np.nanpercentile(mom0, 98)))
    axs[1].imshow(mom0, origin="lower", cmap="RdBu_r", vmin=-mvmax, vmax=mvmax)
    axs[1].set_title(f"mom0 (±{on_kms:.0f} km/s)\nσ={m['sigma_mom0']:.2g}")
    axs[1].set_xticks([]); axs[1].set_yticks([])
    # panel 3: masked mom1
    if np.isfinite(mom1).any():
        v1, v2 = np.nanpercentile(mom1, [5, 95])
        axs[2].imshow(mom1, origin="lower", cmap="RdBu_r", vmin=v1, vmax=v2)
        axs[2].set_title("mom1 (mask mom0>5σ)")
    else:
        axs[2].text(0.5, 0.5, "no mom1\n(no 5σ pixels)", ha="center", va="center",
                     transform=axs[2].transAxes)
        axs[2].set_title("mom1")
    axs[2].set_xticks([]); axs[2].set_yticks([])
    # apply 3x zoom to image panels by setting xlim/ylim around the cutout center
    ny, nx = mom0.shape
    cx, cy = (nx - 1) / 2.0, (ny - 1) / 2.0
    hx, hy = nx / (2.0 * ZOOM), ny / (2.0 * ZOOM)
    for _ax in (axs[0], axs[1], axs[2]):
        _ax.set_xlim(cx - hx, cx + hx)
        _ax.set_ylim(cy - hy, cy + hy)
    # panel 4: spectrum
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
    # panel 5: text summary
    axs[4].axis("off")
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
    axs[4].text(0.05, 0.95, info, transform=axs[4].transAxes, va="top",
                family="monospace", fontsize=9)
    fig.suptitle(f"src {sid} | {line['name']}", fontsize=11)
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


def plot_salt_stack(stacks, sdir, sid):
    if not stacks: return
    fig, axs = plt.subplots(len(stacks), 1, figsize=(10, 3*len(stacks)), sharex=True)
    if len(stacks) == 1: axs = [axs]
    for ax, (species, d) in zip(axs, stacks.items()):
        ax.plot(d["vaxis"], d["spec"], "k-", drawstyle="steps-mid", lw=1)
        ax.axhline(0, color="grey", lw=0.5)
        ax.axhline(d["sigma"], color="red", lw=0.5, ls=":", label=f"1σ={d['sigma']:.2g}")
        ax.axhline(-d["sigma"], color="red", lw=0.5, ls=":")
        ax.axhline(3*d["sigma"], color="red", lw=0.5, ls="--")
        ax.axvline(0, color="blue", lw=0.5, ls="--", label="v_lsr")
        ax.set_ylabel(f"{species} stack (N={d['n_lines']})")
        ax.legend(loc="best", fontsize=8)
    axs[-1].set_xlabel("v - v_lsr (km/s)")
    fig.suptitle(f"Source {sid}: salt stacked spectra")
    fig.tight_layout()
    fig.savefig(sdir / f"salt_stack.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


# ----- Stacking -----
def stack_salt(rows, sdir, vlsr_kms):
    out = {}
    for species in ("NaCl", "KCl"):
        sub = [r for r in rows if r.get("group") == species and r.get("in_band")
               and "peak_Kkms_or_unit" in r]
        if len(sub) < 2:
            continue
        stacks, sigmas = [], []
        v_grid = np.arange(-VWIN_KMS, VWIN_KMS + 0.5, 1.0)
        for r in sub:
            sf = sdir / f"{r['line']}.spec.npz"
            if not sf.exists(): continue
            arr = np.load(sf)
            v_shift = arr["vaxis"] - vlsr_kms
            order = np.argsort(v_shift)
            interp = np.interp(v_grid, v_shift[order], arr["spec"][order],
                                left=np.nan, right=np.nan)
            if not np.isfinite(interp).any():
                continue
            stacks.append(interp); sigmas.append(float(arr["sigma"]))
        if len(stacks) < 2:
            continue
        weights = 1.0 / np.array(sigmas)**2
        stack = np.nansum(np.array(stacks) * weights[:, None], axis=0) / np.nansum(weights)
        sigma_stack = 1.0 / np.sqrt(np.nansum(weights))
        # Reject degenerate stack (all zero or all NaN — means interpolation failed)
        finite_nonzero = np.isfinite(stack) & (stack != 0)
        if finite_nonzero.sum() < 5 or np.nanstd(stack) == 0:
            print(f"  ! {species} stack degenerate (all-zero/NaN); not saving")
            continue
        out[species] = dict(vaxis=v_grid, spec=stack, sigma=sigma_stack, n_lines=len(stacks))
        np.savez(sdir / f"{species}_stack.npz", **out[species])
    return out


# ----- Per-source processing -----
def process_source(src, cubes, lines, outroot, vlsr_kms, distance_kpc, on_kms_default=10.0):
    sid = src["id"]
    sdir = outroot / f"source_{sid:02d}"
    sdir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== Source {sid}: ra={src['ra_deg']:.5f} dec={src['dec_deg']:.5f} "
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
                                     on_kms=70.0 if line["group"]=="RRL" else 15.0)
        hdr_m = cutout.wcs.celestial.to_header()
        hdr_m["BUNIT"] = "km/s"
        fits.PrimaryHDU(m1, hdr_m).writeto(sdir / f"{line['name']}_mom1.fits", overwrite=True)
        hdr0 = hdr_m.copy(); hdr0["BUNIT"] = "Jy/beam km/s"
        fits.PrimaryHDU(m0, hdr0).writeto(sdir / f"{line['name']}_mom0.fits", overwrite=True)
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
    args = ap.parse_args()

    srcdir = ROOT / f"uvdata/{args.proposal}/{args.target}"
    outroot = ROOT / f"analysis_products/{args.target}/{args.proposal}"
    outroot.mkdir(parents=True, exist_ok=True)

    # Filename-target matching to handle multi-source MOUS (e.g., 2022.1.01344.S)
    # ALMA encodes + as p and - as m in filenames.
    target_variants = {args.target,
                       args.target.replace("+", "p").replace("-", "m"),
                       args.target.replace("+", "").replace("-", "")}
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

    sources, cont_sigma, cont_hdr, cont_data, wcs_cont = find_continuum_sources(cont_path, args.peak_sigma)
    print(f"  {len(sources)} sources >= {args.peak_sigma}σ", flush=True)
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
        sdir = outroot / f"source_{int(r['source']):02d}"
        p = sdir / f"{r['line']}_mom0.fits"
        if p.exists():
            top3_paths.append(p)
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
