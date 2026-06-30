"""Publication-quality diagnostics for salt-search detections.

For each (source, family) with a stacked cube under
``analysis_products/<source>/stacks/<family>_stack.fits``, produce:

  1. ``moments/<family>_zoom.png`` — auto-zoomed mom0 + peak + mom1 panel.
     Iteratively lowers the SNR threshold for the moment-1 mask until at
     least *min_pixels* significant pixels are found.  Crops to the bounding
     box of the signal region with a small pad.
  2. ``spectra/<family>_lineid_full.png`` — full-bandwidth raw mean spectrum
     extracted from each contributing native-cube subcube, with NaCl/KCl/H2O
     rest lines overplotted at the source v_LSR.
  3. ``spectra/<family>_stack_velocity.png`` — velocity-stacked spectrum from
     the stack cube, salt lines marked at v=0 (which is the line-stack center).
  4. ``moments/<family>_kepler_score.txt`` — quick metric for whether the
     mom1 map shows a coherent velocity gradient suitable for Keplerian-disk
     mass measurement: max-min spread, smoothness (gradient std), and
     #pixels with S/N>3 in the masked moment0.

Lines are picked from analysis.lines (the all-band list).
"""

import argparse
import json
import logging
import os
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as pl
from astropy import units as u
from astropy.io import fits
from astropy.stats import mad_std
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord
from spectral_cube import SpectralCube
from scipy import ndimage as ndi

from . import paths, lines as line_mod, sources

log = logging.getLogger(__name__)


def _vlsr(name):
    try:
        v = sources.get_vlsr(name).to(u.km / u.s).value
        if np.isfinite(v):
            return v
    except Exception:
        pass
    return 0.0


def _line_color(name):
    if name.startswith("NaCl") or name.startswith("Na37Cl"):
        return "red"
    if name.startswith("KCl") or name.startswith("K37Cl") or name.startswith("KClv") or name.startswith("41KCl"):
        return "blue"
    if name.startswith("H2O"):
        return "darkgreen"
    if name.startswith("SiO"):
        return "orange"
    return "0.4"


def _short_name(name):
    n = re.sub(r"_(\d+)-(\d+)$", r" \1-\2", name)
    n = n.replace("_v", " v=")
    n = n.replace("v=", "v=")
    return n


def auto_zoom_box(mom0, sigma_thresh=3, min_pixels=5, pad_frac=0.3):
    """Find a tight bounding box around significant pixels in *mom0*.

    Iteratively lowers ``sigma_thresh`` until at least *min_pixels* pixels are
    above threshold.  Returns ``(yslice, xslice, used_thresh)`` or None.
    """
    data = np.asarray(mom0).astype(float)
    finite = np.isfinite(data)
    if finite.sum() == 0:
        return None
    noise = mad_std(data[finite])
    if not np.isfinite(noise) or noise <= 0:
        return None
    for thresh in (sigma_thresh, 2.0, 1.5, 1.0, 0.7, 0.5):
        mask = (data > thresh * noise) & finite
        if mask.sum() >= min_pixels:
            ys, xs = np.where(mask)
            ymin, ymax = ys.min(), ys.max()
            xmin, xmax = xs.min(), xs.max()
            dy = max(int((ymax - ymin) * pad_frac), 5)
            dx = max(int((xmax - xmin) * pad_frac), 5)
            ny, nx = data.shape
            return (slice(max(ymin - dy, 0), min(ymax + dy + 1, ny)),
                    slice(max(xmin - dx, 0), min(xmax + dx + 1, nx)),
                    thresh)
    return None


def make_zoom_panel(stack_path, vcen_kms, outpath, family,
                    dv_window=8 * u.km / u.s):
    """Three-panel mom0 / peak / mom1 with auto-zoom."""
    cube = SpectralCube.read(stack_path).with_spectral_unit(
        u.km / u.s, velocity_convention="radio")
    sub = cube.spectral_slab(-dv_window, dv_window)
    sub.allow_huge_operations = True
    mom0 = sub.moment0().value
    peak = sub.max(axis=0).value

    # Iteratively reduce S/N threshold for moment-1 mask
    finite = np.isfinite(peak)
    noise = mad_std(peak[finite]) if finite.any() else 0
    moment1 = None
    used_thresh = None
    for thresh in (5, 3, 2, 1.5, 1.0, 0.5):
        if noise <= 0:
            break
        mask = (peak > thresh * noise) & finite
        if mask.sum() < 3:
            continue
        try:
            m1 = sub.with_mask(mask[None, :, :]).moment1().value
        except Exception:
            continue
        if np.isfinite(m1).sum() < 3:
            continue
        moment1 = m1
        used_thresh = thresh
        break

    box = auto_zoom_box(mom0)
    if box is None:
        ys, xs = slice(None), slice(None)
    else:
        ys, xs, _ = box

    fig = pl.figure(figsize=(15, 5))
    ax1 = fig.add_subplot(1, 3, 1)
    img = mom0[ys, xs]
    norm = simple_norm(img, stretch="asinh", percent=99.5)
    ax1.imshow(img, origin="lower", norm=norm, cmap="viridis")
    ax1.set_title(f"{family} moment 0 (km/s)")
    pl.colorbar(ax1.images[0], ax=ax1, fraction=0.046)

    ax2 = fig.add_subplot(1, 3, 2)
    img2 = peak[ys, xs]
    ax2.imshow(img2, origin="lower",
               norm=simple_norm(img2, stretch="asinh", percent=99.5),
               cmap="inferno")
    ax2.set_title(f"{family} peak (K)")
    pl.colorbar(ax2.images[0], ax=ax2, fraction=0.046)

    ax3 = fig.add_subplot(1, 3, 3)
    if moment1 is not None:
        m1c = moment1[ys, xs]
        finite_m1 = np.isfinite(m1c)
        if finite_m1.sum() > 0:
            v = m1c[finite_m1]
            vmin, vmax = np.nanpercentile(v, [5, 95])
            if vmax - vmin < 1:
                vmin, vmax = vcen_kms - 5, vcen_kms + 5
            ax3.imshow(m1c, origin="lower", cmap="RdBu_r",
                       vmin=vmin, vmax=vmax)
            ax3.set_title(f"{family} mom1 ({used_thresh}σ mask)")
            pl.colorbar(ax3.images[0], ax=ax3, fraction=0.046)
        else:
            ax3.text(0.5, 0.5, "no significant pixels",
                     transform=ax3.transAxes, ha="center", va="center")
    else:
        ax3.text(0.5, 0.5, "mom1 unavailable", transform=ax3.transAxes,
                 ha="center", va="center")
    ax3.set_xticks([]); ax3.set_yticks([])

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    pl.close(fig)
    return outpath, used_thresh


def make_lineid_full(cube_paths, outpath, source_name, vcen_kms,
                     region_radius_arcsec=2.0,
                     peak_coord=None, mode="peak_pixel"):
    """Full-bandwidth spectrum from each native cube, salt lines marked.

    mode='peak_pixel' (default): extract the single-pixel spectrum at
        ``peak_coord`` (preferred) or at ``source_name``'s catalog coord.
        Best for unresolved hot-core lines that get diluted by averaging.
    mode='aperture': mean spectrum over a ``region_radius_arcsec`` circle.
    """
    from astropy.coordinates import SkyCoord
    import regions
    if peak_coord is not None:
        coord = peak_coord
    else:
        src = sources.get_source(source_name)
        coord = SkyCoord(src["RAJ2000"] * u.deg, src["DEJ2000"] * u.deg)
    reg = regions.CircleSkyRegion(coord, radius=region_radius_arcsec * u.arcsec)

    pl.figure(figsize=(16, 4 * max(1, len(cube_paths))))
    n = len(cube_paths)
    if n == 0:
        return None
    for i, cp in enumerate(cube_paths, 1):
        try:
            c = SpectralCube.read(cp)
            if mode == "peak_pixel":
                # 1-pixel spectrum at the peak position.
                wcs_c = c.wcs.celestial
                xp, yp = wcs_c.world_to_pixel(coord)
                xp = int(round(float(xp))); yp = int(round(float(yp)))
                ny, nx = c.shape[1], c.shape[2]
                if 0 <= xp < nx and 0 <= yp < ny:
                    spec = c[:, yp, xp]
                else:
                    spec = None
            else:
                try:
                    sub = c.subcube_from_regions([reg])
                except (ValueError, KeyError, IndexError):
                    ny, nx = c.shape[1], c.shape[2]
                    sub = c[:, max(ny // 2 - 50, 0):ny // 2 + 50,
                            max(nx // 2 - 50, 0):nx // 2 + 50]
                spec = sub.mean(axis=(1, 2))
            if spec is None:
                continue
        except (OSError, ValueError, RuntimeError) as e:
            log.warning("could not extract spectrum from %s: %s", cp, e)
            continue
        freq = spec.spectral_axis.to(u.GHz).value
        ax = pl.subplot(n, 1, i)
        ax.plot(freq, spec.value, color="black", linewidth=0.7)
        ax.set_xlabel("Frequency (GHz, observer)")
        ax.set_ylabel("<S> (Jy/beam or K)")
        f_lo, f_hi = freq.min(), freq.max()

        # Clip y-axis to +20σ / -10σ window to hide deep absorption extremes
        sv = np.asarray(spec.value, dtype=float)
        finite_sv = sv[np.isfinite(sv)]
        if finite_sv.size > 10:
            med = float(np.nanmedian(finite_sv))
            sigma = float(mad_std(finite_sv))
            if np.isfinite(sigma) and sigma > 0:
                ax.set_ylim(med - 10 * sigma, med + 20 * sigma)

        c_kms = 299792.458
        for fam_name, pool in (("NaCl", line_mod.NACL_LINES),
                                ("KCl",  line_mod.KCL_LINES),
                                ("H2O",  line_mod.WATER_LINES),
                                ("ref",  line_mod.REFERENCE_LINES)):
            for lname, lfreq in pool.items():
                f0 = lfreq.to(u.GHz).value * (1 - vcen_kms / c_kms)
                if not (f_lo <= f0 <= f_hi):
                    continue
                col = _line_color(lname)
                ax.axvline(f0, color=col, linestyle=":", linewidth=0.6, alpha=0.7)
                _, ytop = ax.get_ylim()
                ax.text(f0, ytop, _short_name(lname), rotation=90,
                        color=col, fontsize=8, ha="right", va="top",
                        alpha=0.85)
        # Extra ISM + COM/CH3OCHO arrows via the shared style helper.
        try:
            from analysis.lineid_style import apply_labels
            _, ytop = ax.get_ylim()
            apply_labels(ax, f_lo, f_hi, ytop, vsys=vcen_kms,
                          species=("ism",),
                          arrows=("ch3ocho", "ch3oh"), fontsize=8)
        except ImportError:
            pass
        ax.set_xlim(f_lo, f_hi)
        ax.set_title(os.path.basename(cp), fontsize=10)
        ax.tick_params(labelsize=10)

    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    pl.tight_layout()
    pl.savefig(outpath, dpi=300, bbox_inches="tight")
    # Also save a vector PDF companion for paper inclusion.
    base, _ = os.path.splitext(outpath)
    pl.savefig(base + ".pdf", bbox_inches="tight")
    pl.close()
    return outpath


def make_stack_spectrum(stack_path, family, outpath, lines_used):
    """Velocity-stacked spectrum, lines marked at v=0."""
    cube = SpectralCube.read(stack_path)
    cube = cube.with_spectral_unit(u.km / u.s, velocity_convention="radio")
    spec = cube.mean(axis=(1, 2))
    v = spec.spectral_axis.to(u.km / u.s).value

    fig = pl.figure(figsize=(10, 4))
    ax = fig.add_subplot()
    ax.plot(v, spec.value, color="black")
    ax.axvline(0, color="grey", linestyle="--", linewidth=0.6,
               label="line center (stacked)")
    ax.set_xlabel("Velocity offset from rest line (km/s)")
    ax.set_ylabel("<T> (K)")
    ax.set_title(f"{family} stack — lines used: {', '.join(lines_used or [])}")
    ax.legend(loc="upper right", fontsize=8)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    pl.close(fig)
    return outpath


def keplerian_score(mom0_path, peak_path, mom1_path, dist_kpc=None):
    """Quick metric: vrange (max-min in mom1 [km/s]), n_signif_px (>3σ in mom0),
    gradient_std (std of gradient direction) — for ranking sources by
    suitability for Keplerian disk fitting.
    """
    out = {"vrange_kms": np.nan, "n_signif_px": 0, "smoothness": np.nan}
    try:
        mom0 = fits.getdata(mom0_path)
    except Exception:
        return out
    try:
        mom1 = fits.getdata(mom1_path)
    except Exception:
        return out
    finite0 = np.isfinite(mom0)
    if finite0.sum() < 3:
        return out
    noise = mad_std(mom0[finite0])
    if not np.isfinite(noise) or noise <= 0:
        return out
    sig = (mom0 > 3 * noise) & finite0
    out["n_signif_px"] = int(sig.sum())
    if sig.sum() < 5:
        return out
    v = mom1[sig & np.isfinite(mom1)]
    if v.size > 0:
        out["vrange_kms"] = float(np.percentile(v, 95) - np.percentile(v, 5))
    # smoothness: std of central differences inside mask
    gy, gx = np.gradient(np.where(sig, mom1, np.nan))
    g = np.hypot(gy, gx)
    if np.isfinite(g).any():
        out["smoothness"] = float(np.nanstd(g))
    return out


def process_one(source_name, families=("NaCl", "KCl", "H2O")):
    """Generate diagnostics for one survey_2026 source."""
    src_dir = paths.ANALYSIS / source_name
    rep_path = src_dir / "report.json"
    if not rep_path.exists():
        log.warning("no report for %s", source_name)
        return None
    with open(rep_path) as fh:
        rep = json.load(fh)
    vcen = rep.get("vcen_kms", 0.0) or 0.0
    cube_paths = rep.get("cube_paths", []) or []

    out = {"source": source_name, "diagnostics": {}}

    # full-bandwidth line-ID spectrum (one PNG, all cubes stacked vertically)
    if cube_paths:
        lineid_png = src_dir / "spectra" / "lineid_full.png"
        try:
            # make_lineid_full saves both .png and .pdf for paper inclusion.
            make_lineid_full(cube_paths, str(lineid_png), source_name, vcen)
            out["lineid_full"] = str(lineid_png)
            out["lineid_full_pdf"] = str(lineid_png).replace(".png", ".pdf")
        except Exception as e:
            log.warning("lineid_full failed: %s", e)

    for fam in families:
        stack = src_dir / "stacks" / f"{fam}_stack.fits"
        if not stack.exists():
            continue
        # zoom panel
        zoom_png = src_dir / "moments" / f"{fam}_zoom.png"
        try:
            make_zoom_panel(str(stack), vcen, str(zoom_png), fam)
        except Exception as e:
            log.warning("%s %s zoom failed: %s", source_name, fam, e)
        # stack-spectrum with lines used
        stk_png = src_dir / "spectra" / f"{fam}_stack_velocity.png"
        try:
            lines_used = (rep.get("stacks", {}).get(fam, {}) or {}).get("lines_used")
            make_stack_spectrum(str(stack), fam, str(stk_png), lines_used)
        except Exception as e:
            log.warning("%s %s stack-spec failed: %s", source_name, fam, e)
        out["diagnostics"][fam] = {"zoom": str(zoom_png), "stack_spectrum": str(stk_png)}
    return out


def kepler_assess(families=("NaCl", "KCl", "H2O"),
                  out_csv=None,
                  min_peak_thresh=3.0,
                  min_signif_px=20,
                  vrange_min_kms=3.0,
                  vrange_max_kms=25.0):
    """Walk all stacked sources; rank by velocity-gradient strength.

    Reject pure-noise mom1 maps. Require:
      - peak mask >= min_peak_thresh sigma (no fallback below)
      - at least min_signif_px px in mom0 > 5σ
      - vrange in [vrange_min, vrange_max] km/s
      - vrange < 0.9 × full slab width (else sign of pure-noise mom1)
    Score: peak SNR-weighted vrange; coherent gradient bonus.
    """
    from astropy.table import Table
    rows = []
    for d in sorted(paths.ANALYSIS.glob("*/")):
        rep_path = d / "report.json"
        if not rep_path.exists():
            continue
        with open(rep_path) as fh:
            rep = json.load(fh)
        for fam in families:
            stack = d / "stacks" / f"{fam}_stack.fits"
            if not stack.exists():
                continue
            try:
                cube = SpectralCube.read(str(stack)).with_spectral_unit(
                    u.km / u.s, velocity_convention="radio")
                slab_halfwidth = 8.0  # km/s
                sub = cube.spectral_slab(-slab_halfwidth * u.km / u.s,
                                         slab_halfwidth * u.km / u.s)
                slab_full_width = 2 * slab_halfwidth
                mom0 = sub.moment0().value
                peak = sub.max(axis=0).value

                fin_peak = np.isfinite(peak)
                if fin_peak.sum() < 10:
                    continue
                peak_noise = mad_std(peak[fin_peak])
                if not np.isfinite(peak_noise) or peak_noise <= 0:
                    continue
                peak_snr_max = float(np.nanmax(peak[fin_peak]) / peak_noise)
                if peak_snr_max < 5:
                    # nothing detected at >5σ peak — skip
                    continue

                # Hard threshold (no <2σ fallback)
                mask = peak > min_peak_thresh * peak_noise
                if mask.sum() < min_signif_px:
                    continue

                m1 = sub.with_mask(mask[None, :, :]).moment1().value
                finite = np.isfinite(m1) & np.isfinite(mom0) & mask
                if finite.sum() < min_signif_px:
                    continue
                v = m1[finite]
                vrange = float(np.percentile(v, 95) - np.percentile(v, 5))

                # Reject pure-noise mom1: vrange spanning > 90% slab
                if vrange >= 0.9 * slab_full_width:
                    continue
                # Reject too-narrow (no rotation) or too-wide (contamination)
                if not (vrange_min_kms <= vrange <= vrange_max_kms):
                    continue

                # mom0 SNR
                fin0 = np.isfinite(mom0)
                mom0_noise = mad_std(mom0[fin0])
                if not np.isfinite(mom0_noise) or mom0_noise <= 0:
                    continue
                nsig = int(((mom0 > 5 * mom0_noise) & fin0).sum())
                if nsig < min_signif_px:
                    continue
                mom0_snr_max = float(np.nanmax(mom0[fin0]) / mom0_noise)

                # Coherent-gradient check: dipole orientation consistency.
                # Split mask into low-velocity-half vs high-velocity-half by
                # median; require they be spatially separated (centroids apart).
                vmed = float(np.median(v))
                lo_mask = (m1 < vmed) & finite
                hi_mask = (m1 >= vmed) & finite
                if lo_mask.sum() < 3 or hi_mask.sum() < 3:
                    continue
                ys, xs = np.indices(m1.shape)
                lo_cx, lo_cy = xs[lo_mask].mean(), ys[lo_mask].mean()
                hi_cx, hi_cy = xs[hi_mask].mean(), ys[hi_mask].mean()
                centroid_sep_px = float(np.hypot(hi_cx - lo_cx, hi_cy - lo_cy))
                # require centroids separated by ≥ ~1.5 px
                if centroid_sep_px < 1.5:
                    continue

                # gradient smoothness
                gy, gx = np.gradient(np.where(finite, m1, np.nan))
                gmag = np.hypot(gy, gx)
                grad_med = float(np.nanmedian(gmag))

                # Score: rewards peak SNR × vrange × coherent gradient.
                # log10(nsig) bounded so it doesn't dominate.
                score = (vrange
                         * np.log10(max(nsig, 1))
                         * np.log10(max(peak_snr_max, 1.01))
                         * np.log10(max(centroid_sep_px, 1.01)))

                rows.append({
                    "source": d.name, "family": fam,
                    "n_signif_px": nsig,
                    "peak_snr_max": peak_snr_max,
                    "mom0_snr_max": mom0_snr_max,
                    "vrange_kms": vrange,
                    "centroid_sep_px": centroid_sep_px,
                    "grad_median": grad_med,
                    "kepler_score": float(score),
                })
            except Exception as e:
                log.warning("kepler eval failed for %s/%s: %s", d.name, fam, e)
    if not rows:
        print("no rows")
        return
    t = Table(rows=rows)
    t.sort("kepler_score")
    t.reverse()
    if out_csv:
        t.write(out_csv, format="csv", overwrite=True)
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sources", nargs="*",
                    help="explicit source names; default = every survey_2026 dir with stacks")
    ap.add_argument("--kepler", action="store_true",
                    help="compute Keplerian-suitability scoring table only")
    ap.add_argument("--out-csv", default=str(paths.DATA / "kepler_candidates.csv"))
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    if args.kepler:
        t = kepler_assess(out_csv=args.out_csv)
        if t is not None:
            print(t[:30])
        return

    if args.sources:
        targets = args.sources
    else:
        targets = sorted(d.name for d in paths.ANALYSIS.glob("*/")
                         if (d / "report.json").exists()
                         and any(d.glob("stacks/*_stack.fits")))
    print(f"processing {len(targets)} sources")
    for t in targets:
        try:
            process_one(t)
            print(f"  {t}: done")
        except Exception as e:
            log.exception("%s failed: %s", t, e)


if __name__ == "__main__":
    main()
