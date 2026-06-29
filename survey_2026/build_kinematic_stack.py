"""Kinematic-guided per-cube stack for any target.

Given a `target`, `proposal`, and a `guide_line` rest frequency, build a
mom1 of the guide line and use it as the velocity surface to shift every
other cube's spectrum onto the same rest-velocity grid. Output one combined
'aligned' spectrum spanning the full observed frequency range, plus a labeled
line-id plot with NaCl/KCl/H2O/RRL/COM transitions.

Usage:
    python build_kinematic_stack.py --target MonR2-IRS2 --proposal 2018.1.00446.S \
        --guide-line H26alpha --vlsr 10.0

For every L4_d2 source we want at least one of these stacks. The brightest
disk-associated line is the right guide because the velocity field traces
the disk gradient that shifts every other line by the same amount.
"""
import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import mad_std
from astropy.wcs import WCS
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
C_KMS = 299792.458

# Common line-id list: rest GHz, label
LINELIST = [
    (217.10498, r"SiO 5-4"),
    (217.97715, r"NaCl$_{v=0}$ 17-16"),
    (218.22219, r"H$_2$CO 3$_{0,3}$-2$_{0,2}$"),
    (219.56035, r"C$^{18}$O 2-1"),
    (219.61300, r"NaCl$_{v=1}$ 17-16"),
    (220.39870, r"$^{13}$CO 2-1"),
    (221.96527, r"H30$\alpha$"),
    (230.53800, r"CO 2-1"),
    (230.77580, r"NaCl$_{v=2}$ 18-17"),
    (231.06099, r"OCS 19-18"),
    (231.90050, r"H30$\alpha$"),
    (232.18712, r"$^{13}$CS 5-4"),
    (232.68673, r"H$_2$O 5$_{5,0}$-6$_{4,3}$"),
    (232.93657, r"H$_2$O$_{v_2=1}$ 3$_{1,3}$-2$_{2,0}$"),
    (234.14050, r"H$_2$O$_{v_2=1}$ 4$_{1,4}$-3$_{2,1}$"),
    (234.25090, r"NaCl$_{v=0}$ 18-17"),
    (234.41450, r"H$_2$O$_{v_2=1}$ 9$_{1,8}$-8$_{4,5}$"),
]

# Common RRLs
RRLS = {
    "H26alpha": 353.6228,
    "H27alpha": 316.4156,
    "H28alpha": 283.5471,
    "H29alpha": 256.3022,
    "H30alpha": 231.9009,
    "H31alpha": 210.5018,
    "H32alpha": 192.0193,
    "H39beta": 308.1064,
    "H40beta": 286.5181,
}


def jy_per_beam_to_K(nu_GHz, bmaj_arcsec, bmin_arcsec):
    return 1.222e6 / (nu_GHz ** 2 * bmaj_arcsec * bmin_arcsec)


def cube_freq_range_GHz(hdr):
    nu0 = float(hdr["CRVAL3"])
    dnu = float(hdr["CDELT3"])
    crp = float(hdr["CRPIX3"])
    n = int(hdr["NAXIS3"])
    f = (nu0 + (np.arange(n) + 1 - crp) * dnu) / 1e9
    return f, dnu / 1e9


def find_cube_with_line(cubes, rest_GHz, vlsr_kms):
    f_obs = rest_GHz * (1.0 - vlsr_kms / C_KMS)
    for c in cubes:
        h = fits.getheader(c)
        f, _ = cube_freq_range_GHz(h)
        if f.min() + 0.05 <= f_obs <= f.max() - 0.05:
            return c
    return None


def mom1_at_pixel(cube_path, src_pix, rest_GHz, vlsr_kms, half_kms=30.0,
                    thresh_sigma=3.0):
    """Build channel-of-peak velocity field within ~5 beam around src; return
    mom1 (km/s) at the source pixel, or None if S/N too low."""
    with fits.open(cube_path, memmap=True) as hdul:
        h = hdul[0].header
        data = hdul[0].data
    if data.ndim == 4:
        data = data[0]
    f, _ = cube_freq_range_GHz(h)
    v = (1.0 - f / rest_GHz) * C_KMS
    on = np.abs(v - vlsr_kms) < half_kms
    if on.sum() < 3:
        return None
    sx, sy = src_pix
    nx, ny = data.shape[2], data.shape[1]
    half_pix = 30
    x0 = max(0, sx - half_pix); x1 = min(nx, sx + half_pix + 1)
    y0 = max(0, sy - half_pix); y1 = min(ny, sy + half_pix + 1)
    sub = data[on][:, y0:y1, x0:x1]
    sub = np.asarray(sub, dtype=np.float32)
    sigma = float(mad_std(sub, ignore_nan=True))
    weight = np.clip(sub - thresh_sigma * sigma, 0, None)
    if weight.sum() <= 0:
        return None
    v_sub = v[on][:, None, None]
    mom1 = np.nansum(weight * v_sub, axis=0) / np.nansum(weight, axis=0)
    cy = sy - y0; cx = sx - x0
    val = mom1[cy, cx]
    if not np.isfinite(val):
        return None
    return float(val)


def extract_aligned_spectrum(cube_path, src_pix, mom1_kms):
    """Extract source-pixel spectrum (in K), shifted so the rest velocity
    frame is centered at v=0 using mom1_kms as the per-cube velocity offset."""
    with fits.open(cube_path, memmap=True) as hdul:
        h = hdul[0].header
        data = hdul[0].data
    if data.ndim == 4:
        data = data[0]
    sx, sy = src_pix
    if not (0 <= sx < data.shape[2] and 0 <= sy < data.shape[1]):
        return None
    spec_jy = np.asarray(data[:, sy, sx], dtype=np.float32)
    f, _ = cube_freq_range_GHz(h)  # GHz
    bmaj = float(h.get("BMAJ", 0)) * 3600
    bmin = float(h.get("BMIN", 0)) * 3600
    K = jy_per_beam_to_K(float(np.median(f)), bmaj, bmin)
    spec_K = spec_jy * K
    # Shift observed frequency by mom1_kms back to rest frame:
    # f_rest = f_obs / (1 - v_obs/c). f_obs here is f from the header.
    f_rest_GHz = f / (1.0 - mom1_kms / C_KMS)
    return f_rest_GHz, spec_K


def overlay_lineids(ax, freq_GHz_data, flux_data):
    try:
        from lineid_plot import plot_line_ids
    except ImportError:
        return
    fmin, fmax = float(freq_GHz_data.min()), float(freq_GHz_data.max())
    keep_freqs = []; keep_labels = []
    for rest, lab in LINELIST:
        if fmin <= rest <= fmax:
            keep_freqs.append(rest); keep_labels.append(lab)
    for n, rest in RRLS.items():
        if fmin <= rest <= fmax:
            keep_freqs.append(rest); keep_labels.append(n)
    if not keep_freqs:
        return
    plot_line_ids(np.asarray(freq_GHz_data, dtype=float),
                   np.asarray(flux_data, dtype=float),
                   np.array(keep_freqs, dtype=float),
                   keep_labels, ax=ax,
                   label1_size=8, extend=False, max_iter=200)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--proposal", required=True)
    ap.add_argument("--guide-line", required=True,
                     help="One of 'H26alpha', 'H30alpha', 'NaCl_v0_J18-17', "
                          "or 'auto' (= brightest line in measurements).")
    ap.add_argument("--vlsr", type=float, required=True)
    ap.add_argument("--source-id", type=int, default=None,
                     help="If omitted, uses brightest mm continuum source.")
    args = ap.parse_args()

    prop_dir = ROOT / "analysis_products" / args.target / args.proposal
    if not prop_dir.is_dir():
        raise SystemExit(f"no analysis_products dir: {prop_dir}")
    cont = pd.read_csv(prop_dir / "continuum_sources.csv")
    if args.source_id is None:
        bid = int(cont.loc[cont["peak_Jybeam"].idxmax(), "id"])
    else:
        bid = args.source_id
    ra = float(cont.loc[bid - 1, "ra_deg"])
    dec = float(cont.loc[bid - 1, "dec_deg"])
    coord = SkyCoord(ra, dec, unit="deg")
    print(f"target={args.target} src #{bid} at {ra:.5f},{dec:.5f}, vlsr={args.vlsr}")

    cubes = sorted((ROOT / "uvdata" / args.proposal / args.target)
                     .glob("*cube*pbcor.fits"))
    if not cubes:
        raise SystemExit(f"no cubes for {args.target}/{args.proposal}")

    # Resolve guide line rest freq
    guide_rest = None
    if args.guide_line in RRLS:
        guide_rest = RRLS[args.guide_line]
    if guide_rest is None:
        rest_lookup = {l[1].split()[0]: l[0] for l in LINELIST}
        rest_lookup.update({n: f for n, f in RRLS.items()})
        guide_rest = rest_lookup.get(args.guide_line)
    if guide_rest is None:
        # Last resort: pull rest_GHz from line_measurements.csv
        meas = prop_dir / "line_measurements.csv"
        if meas.exists():
            mdf = pd.read_csv(meas)
            row = mdf[mdf["line"] == args.guide_line]
            if not row.empty and np.isfinite(row.iloc[0].get("rest_GHz", np.nan)):
                guide_rest = float(row.iloc[0]["rest_GHz"])
    if guide_rest is None:
        raise SystemExit(f"don't know rest freq for {args.guide_line}; "
                          f"add to RRLS / LINELIST or to line_measurements.csv")
    guide_cube = find_cube_with_line(cubes, guide_rest, args.vlsr)
    if guide_cube is None:
        raise SystemExit(f"no cube contains guide line {args.guide_line} "
                          f"({guide_rest} GHz)")
    print(f"guide line {args.guide_line} ({guide_rest} GHz) in {guide_cube.name}")

    # Source pixel in the guide cube
    h_g = fits.getheader(guide_cube)
    wcs_g = WCS(h_g).celestial
    xp, yp = wcs_g.world_to_pixel(coord)
    src_pix_g = (int(round(float(xp))), int(round(float(yp))))
    mom1_g = mom1_at_pixel(guide_cube, src_pix_g, guide_rest, args.vlsr)
    if mom1_g is None:
        print("guide line too weak to build mom1; falling back to vlsr")
        mom1_g = args.vlsr
    print(f"guide mom1 at source pixel: {mom1_g:.2f} km/s")

    # Extract aligned spectra from every cube; keep per-cube so we can panel
    panels = []  # list of dicts {mous, spw, freq_rest, T_K}
    for cube in cubes:
        h = fits.getheader(cube)
        wcs = WCS(h).celestial
        xp2, yp2 = wcs.world_to_pixel(coord)
        ix = int(round(float(xp2))); iy = int(round(float(yp2)))
        if not (0 <= ix < h["NAXIS1"] and 0 <= iy < h["NAXIS2"]):
            print(f"  skip {cube.name}: source off-FOV")
            continue
        out = extract_aligned_spectrum(cube, (ix, iy), mom1_g)
        if out is None:
            continue
        f_rest, T = out
        # MOUS + SPW tag for panel labeling
        stem = cube.stem
        spw = next((t for t in stem.split(".") if t.startswith("spw")),
                    "spw??")
        mous = "obs"
        if stem.startswith("member"):
            for tok in reversed(stem.split(".")[0].split("_")):
                if tok.startswith("X") and len(tok) > 1:
                    mous = tok; break
        panels.append(dict(mous=mous, spw=spw, freq=f_rest, T=T,
                              name=cube.name))
        print(f"  {cube.name}: nu_rest={f_rest.min():.3f}-{f_rest.max():.3f} GHz "
              f"T_max={np.nanmax(T):.1f} K")

    if not panels:
        raise SystemExit("no spectra extracted")

    out_dir = ROOT / "analysis_products" / args.target / args.proposal / "kinematic_stack"
    out_dir.mkdir(exist_ok=True)
    # Save concatenated + per-panel data for downstream use
    f_arr = np.concatenate([p["freq"] for p in panels])
    T_arr = np.concatenate([p["T"] for p in panels])
    order = np.argsort(f_arr)
    f_arr = f_arr[order]; T_arr = T_arr[order]
    np.savez(out_dir / f"aligned_by_{args.guide_line}.npz",
              freq_GHz=f_arr, T_K=T_arr,
              guide_line=args.guide_line, guide_rest_GHz=guide_rest,
              mom1_kms=mom1_g)
    print(f"wrote {out_dir}/aligned_by_{args.guide_line}.npz")

    # Sort panels by MOUS then SPW number for clean ordering
    def _spw_num(s):
        import re as _re
        m = _re.search(r"(\d+)", s)
        return int(m.group(1)) if m else 99
    panels.sort(key=lambda p: (p["mous"], _spw_num(p["spw"])))

    # Paginate at most ROWS_PER_PAGE rows per PNG
    ROWS_PER_PAGE = 8
    n_pages = (len(panels) + ROWS_PER_PAGE - 1) // ROWS_PER_PAGE
    for pi in range(n_pages):
        page = panels[pi * ROWS_PER_PAGE:(pi + 1) * ROWS_PER_PAGE]
        n = len(page)
        fig, axes = plt.subplots(n, 1, figsize=(18, 2.5 * n + 0.6))
        if n == 1:
            axes = [axes]
        for ax, p in zip(axes, page):
            freq = p["freq"]; T = p["T"]
            ax.plot(freq, T, "k-", lw=0.6)
            ax.set_ylabel(f"{p['mous']}\n{p['spw']}\nT (K)", fontsize=7)
            ax.tick_params(labelsize=8)
            ax.set_xlim(float(freq.min()), float(freq.max()))
            finite = T[np.isfinite(T)]
            if finite.size < 5:
                continue
            ymin = float(np.nanpercentile(finite, 1))
            ymax = float(np.nanpercentile(finite, 99)) * 1.15 + 1
            if not (np.isfinite(ymin) and np.isfinite(ymax)) or ymax <= ymin:
                continue
            ax.set_ylim(ymin, ymax)
            overlay_lineids(ax, freq, T)
        axes[-1].set_xlabel("rest frequency (GHz, aligned)")
        suffix = "" if n_pages == 1 else f"_p{pi+1}"
        fig.suptitle(
            f"{args.target} src{bid:02d} aligned by {args.guide_line}  "
            f"(mom1={mom1_g:.1f} km/s, vlsr={args.vlsr:+.1f}"
            f"{f', page {pi+1}/{n_pages}' if n_pages > 1 else ''})",
            fontsize=10)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        out_png = out_dir / f"aligned_by_{args.guide_line}{suffix}.png"
        fig.savefig(out_png, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {out_png} ({n} panels)")


if __name__ == "__main__":
    main()
