"""For every analyzed source, write an N-panel SPW plot of the
source-region-averaged spectrum with lineid annotations from the 2019/2023
papers' disk_lines + absorbers dictionaries.

Source region preference order:
  1) sources.reg ellipse/circle around the source from line_pipeline
  2) bright continuum peak pixel + half-beam radius
  3) just the central pixel

For each analyzed (target, proposal, source_id) write:
  analysis_products/<target>/<proposal>/spectrum_panels_src<NN>.png

Lineid catalog source: /orange/adamginsburg/salt/Orion_ALMA_2016.1.00165.S/analysis/lines.py
  (disk_lines + absorbers = 161 entries spanning B3/B6/B7)
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

sys.path.insert(0, "/orange/adamginsburg/salt/Orion_ALMA_2016.1.00165.S/analysis")
import lines as _lines

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"

DISK_LINES = {**_lines.disk_lines, **_lines.absorbers}
LINE_REST_GHZ = {n: float(u.Quantity(f).to(u.GHz).value) for n, f in DISK_LINES.items()}


def jy_per_beam_to_K(nu_GHz, bmaj_arcsec, bmin_arcsec):
    return 1.222e6 / (nu_GHz ** 2 * bmaj_arcsec * bmin_arcsec)


def latex_label(name):
    """Best-effort LaTeX-ification of disk_line names. Falls back to as-is."""
    # Common substitutions
    s = name
    s = s.replace("v=", "$_{v=").replace("_", "}\\,") if "v=" in name else s
    if "$_{v=" in s and not s.endswith("}"):
        # Reclose any open subscript
        s = s + "}"
    return s


def per_source_pixel(meas_csv, sid, target, proposal):
    """Return (ra, dec) of the SOURCE position from continuum_sources.csv."""
    cont = pd.read_csv(ANALYSIS / target / proposal / "continuum_sources.csv")
    row = cont[cont["id"] == sid]
    if row.empty:
        return None
    return float(row.iloc[0]["ra_deg"]), float(row.iloc[0]["dec_deg"])


def avg_spectrum_over_beam(cube_path, coord, n_half_pixels=2):
    """Average spectrum over a (2*n_half_pixels+1)^2 box around the source.

    For unresolved disks the source FWHM is ~1 beam; we mean ~5x5 pixels =
    a few beams. Output: (freq_GHz, T_K)."""
    with fits.open(cube_path, memmap=True) as hdul:
        h = hdul[0].header
        data = hdul[0].data
    if data.ndim == 4:
        data = data[0]
    wcs = WCS(h).celestial
    xp, yp = wcs.world_to_pixel(coord)
    ix = int(round(float(xp))); iy = int(round(float(yp)))
    nx = data.shape[2]; ny = data.shape[1]
    if not (0 <= ix < nx and 0 <= iy < ny):
        return None
    y0 = max(0, iy - n_half_pixels); y1 = min(ny, iy + n_half_pixels + 1)
    x0 = max(0, ix - n_half_pixels); x1 = min(nx, ix + n_half_pixels + 1)
    cutout = data[:, y0:y1, x0:x1]
    cutout = np.asarray(cutout, dtype=np.float32)
    spec_jybeam = np.nanmean(cutout, axis=(1, 2))
    nu0 = float(h.get("CRVAL3"))
    dnu = float(h.get("CDELT3"))
    crp = float(h.get("CRPIX3"))
    nchan = data.shape[0]
    freq_Hz = nu0 + (np.arange(nchan) + 1 - crp) * dnu
    bmaj = float(h.get("BMAJ", 0)) * 3600
    bmin = float(h.get("BMIN", 0)) * 3600
    nu_GHz = float(np.median(freq_Hz)) / 1e9
    K = jy_per_beam_to_K(nu_GHz, bmaj, bmin)
    return freq_Hz / 1e9, spec_jybeam * K


def overlay_lineids_panel(ax, fmin, fmax, ymax):
    keep = [(n, f) for n, f in LINE_REST_GHZ.items()
            if fmin <= f <= fmax]
    if not keep:
        return
    try:
        from lineid_plot import plot_line_ids
    except ImportError:
        for _, f in keep:
            ax.axvline(f, color="C1", lw=0.5, alpha=0.5)
        return
    keep.sort(key=lambda kv: kv[1])
    names = [n for n, _ in keep]
    freqs = np.array([f for _, f in keep], dtype=float)
    # Build a dummy 'data' for label placement matching y range
    dummy_wave = np.linspace(fmin, fmax, 200)
    dummy_flux = np.full(200, ymax)
    plot_line_ids(dummy_wave, dummy_flux, freqs, names, ax=ax,
                   label1_size=6, extend=False, max_iter=50)


def plot_source(target, proposal, src_id):
    """Make the N-panel SPW spectrum plot for one source."""
    target_dir = ANALYSIS / target / proposal
    coord = per_source_pixel(None, src_id, target, proposal)
    if coord is None:
        print(f"  no continuum source row for src {src_id}")
        return
    ra, dec = coord
    skycoord = SkyCoord(ra, dec, unit="deg")

    uvdir = ROOT / "uvdata" / proposal / target
    cubes = sorted(uvdir.glob("*sci*cube*pbcor.fits"))
    if not cubes:
        print(f"  no cubes for {target}/{proposal}")
        return

    panels = []
    for cube in cubes:
        out = avg_spectrum_over_beam(cube, skycoord)
        if out is None:
            continue
        spw = next((t for t in cube.stem.split(".") if t.startswith("spw")),
                    "spw??")
        panels.append((spw, *out))
    if not panels:
        print(f"  source off-FOV in all cubes ({target} src{src_id:02d})")
        return

    n = len(panels)
    fig, axes = plt.subplots(n, 1, figsize=(12, 1.6 * n + 0.6))
    if n == 1:
        axes = [axes]
    for ax, (spw, freq, T) in zip(axes, panels):
        ax.plot(freq, T, "k-", lw=0.6)
        ax.set_ylabel(f"{spw}\nT (K)", fontsize=8)
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
        overlay_lineids_panel(ax, float(freq.min()), float(freq.max()), ymax)
    axes[-1].set_xlabel("rest frequency [observed] (GHz)")
    fig.suptitle(f"{target} src{src_id:02d}  (proposal {proposal})", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out_png = target_dir / f"spectrum_panels_src{src_id:02d}.png"
    fig.savefig(out_png, dpi=110)
    plt.close(fig)
    print(f"  wrote {out_png}")


def main():
    # Audit gives us the canonical (target, proposal, src_id) tuples
    audit = pd.read_csv(ROOT / "data/evidence_audit.csv")
    if audit.empty:
        print("nothing in audit")
        return
    for _, r in audit.iterrows():
        print(f"==> {r['target']} {r['proposal']} src{int(r['src_id']):02d}")
        plot_source(str(r["target"]), str(r["proposal"]), int(r["src_id"]))


if __name__ == "__main__":
    main()
