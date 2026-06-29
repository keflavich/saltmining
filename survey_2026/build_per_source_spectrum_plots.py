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
import re
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

# Restricted subset used for the per-panel overlay: salts + the strongest
# COMs + common shock/HII tracers in the band. Keeps the panels readable
# (lineid_plot fails to lay out 60+ labels in a 12" panel).
_KEEP_PATTERNS = ("NaCl", "KCl", "Na37Cl", "K37Cl", "H2O",
                  "SiO", "29SiO", "SiS", "SO_", "SO2",
                  "H30alpha", "H29alpha", "H26alpha", "H37beta", "H42alpha",
                  "CO_2-1", "13CO", "C18O", "CO-18", "12CO",
                  "CH3OH_", "CH3CN", "OCS", "HC3N", "H2CO303-202",
                  "CH3OCHO", "HC(O)NH2")


def _line_is_important(name):
    return any(p in name for p in _KEEP_PATTERNS)


def jy_per_beam_to_K(nu_GHz, bmaj_arcsec, bmin_arcsec):
    if not (nu_GHz > 0 and bmaj_arcsec > 0 and bmin_arcsec > 0):
        return np.nan
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
    """Simple vertical-line + small rotated text overlay.

    Filters to the salts + strong COMs + RRLs + common shock/HII tracers
    (about 30 lines vs the full 161-line catalog) so the panel stays
    readable. Labels are drawn rotated 90 deg above the spectrum range."""
    keep = [(n, f) for n, f in LINE_REST_GHZ.items()
            if fmin <= f <= fmax and _line_is_important(n)]
    if not keep:
        return
    keep.sort(key=lambda kv: kv[1])
    for n, f in keep:
        ax.axvline(f, color="C1", lw=0.4, alpha=0.55)
        ax.text(f, ymax * 0.97, n, rotation=90, ha="center", va="top",
                 fontsize=5, color="C1", alpha=0.9, clip_on=True)


def mous_tag(cube_path):
    """Short MOUS id from filename: member.uid___A001_X1465_X2fd3.* -> X2fd3."""
    stem = cube_path.stem
    # Patterns: member.uid___A001_X<mous_top>_X<mous_id>.<targetname>_sci.spwNN...
    parts = stem.split(".")
    if parts and parts[0].startswith("member"):
        # e.g. member.uid___A001_X1465_X2fd3 -> grab last X<id>
        head = parts[0]
        sub = head.split("_")
        for tok in reversed(sub):
            if tok.startswith("X") and len(tok) > 1:
                return tok
    return "obs"


def plot_source(target, proposal, src_id, rows_per_page=8):
    """Make N-panel SPW spectrum plot(s) for one source. One panel per
    (MOUS, SPW) cube; paginate to keep at most rows_per_page rows per PNG."""
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
        mous = mous_tag(cube)
        # Sort key: SPW number for clean ordering
        spw_num = int(re.sub(r"\D", "", spw)) if any(c.isdigit() for c in spw) else 99
        panels.append((mous, spw, spw_num, *out))
    if not panels:
        print(f"  source off-FOV in all cubes ({target} src{src_id:02d})")
        return

    # Sort: by MOUS then SPW so observations stay grouped
    panels.sort(key=lambda p: (p[0], p[2]))

    # Paginate
    n_pages = (len(panels) + rows_per_page - 1) // rows_per_page
    for pi in range(n_pages):
        page = panels[pi * rows_per_page:(pi + 1) * rows_per_page]
        n = len(page)
        fig, axes = plt.subplots(n, 1, figsize=(18, 2.5 * n + 0.6))
        if n == 1:
            axes = [axes]
        for ax, (mous, spw, _, freq, T) in zip(axes, page):
            ax.plot(freq, T, "k-", lw=0.6)
            ax.set_ylabel(f"{mous}\n{spw}\nT (K)", fontsize=7)
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
        axes[-1].set_xlabel("observed frequency (GHz)")
        suffix = "" if n_pages == 1 else f"_p{pi+1}"
        fig.suptitle(
            f"{target} src{src_id:02d}  (proposal {proposal}"
            f"{f', page {pi+1}/{n_pages}' if n_pages > 1 else ''})", fontsize=10)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        out_png = target_dir / f"spectrum_panels_src{src_id:02d}{suffix}.png"
        fig.savefig(out_png, dpi=110)
        plt.close(fig)
        print(f"  wrote {out_png} ({n} panels)")


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
