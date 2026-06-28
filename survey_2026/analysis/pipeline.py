"""Top-level driver: ``run_one_source(name, vcen=…)`` runs the full salt
search pipeline on a single source.

Design follows `miriam/notebooks/Sanhueza18162.ipynb`:

  1. Look up source coords / distance / vlsr.
  2. Discover cubes (survey_2026 + DIHCA2).
  3. For each line family (NaCl, KCl, H2O):
       a. open cubes that cover any line in the family;
       b. cut to a small region around the source (``region_radius``);
       c. continuum-subtract with `cube.percentile(25)`;
       d. convolve all to a common beam;
       e. stack via `spectral_cube.analysis_utilities.stack_cube`;
       f. write the stack and a 3-panel quicklook PNG.
  4. Build moment / RGB / PV figures.

Output structure (under ``analysis_products/<source>/``)::

    <source>/
        cubes/         # convolved, region-cut, contsub'd cubes
        stacks/        # NaCl / KCl / H2O stacked cubes (FITS)
        moments/       # moment-0 / RGB / mom1 PNGs
        spectra/       # region-averaged spectra
        pv/            # PV plots
        report.json    # source params and product paths
"""

import json
import logging
from pathlib import Path

import numpy as np
import regions
from astropy import units as u
from astropy.coordinates import SkyCoord

from . import paths, sources, cubes as cube_mod, lines as line_mod
from . import stacking, moments as mom_mod, pv as pv_mod


log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(name)s: %(message)s")


def _detection_stats(stack, vcen, on_window_kms=5, off_min_kms=15,
                      source_coord=None):
    """Peak / noise / SNR for a stacked cube around the source.

    On-source velocity slab: |v - vcen| < on_window_kms.
    Off-line noise: mad_std of pixels with |v - vcen| > off_min_kms.

    Two SNRs are reported:
      snr_global: brightest pixel anywhere in the on-velocity slab (legacy).
      snr_center: spectrum at the source pixel (when source_coord supplied),
                    taking the on-window peak.

    snr_center is the right physical quantity for disk emission centered at
    the YSO; snr_global flags any feature in the field of view (often noise
    spikes at the edge or off-source emission).
    """
    from astropy.stats import mad_std
    from spectral_cube import SpectralCube
    if isinstance(stack, str):
        stack = SpectralCube.read(stack)
    cube = stack.with_spectral_unit(u.km / u.s, velocity_convention="radio")
    v = cube.spectral_axis
    on_mask = np.abs(v - vcen) < on_window_kms * u.km / u.s
    off_mask = np.abs(v - vcen) > off_min_kms * u.km / u.s
    if on_mask.sum() == 0 or off_mask.sum() == 0:
        return {"peak_K": np.nan, "noise_K": np.nan, "snr": np.nan,
                "snr_global": np.nan, "snr_center": np.nan,
                "peak_K_center": np.nan, "peak_v_kms": np.nan}
    on_data = cube.filled_data[on_mask, :, :].value
    off_data = cube.filled_data[off_mask, :, :].value
    peak = float(np.nanmax(on_data))
    noise = float(mad_std(off_data, ignore_nan=True))
    spec = np.nanmax(np.nanmax(on_data, axis=2), axis=1)
    iv = int(np.nanargmax(spec))
    v_at_peak = float(v[on_mask][iv].value)
    snr_global = peak / noise if noise > 0 else np.nan

    peak_center = np.nan; snr_center = np.nan
    if source_coord is not None:
        wcs = cube.wcs.celestial
        xp, yp = wcs.world_to_pixel(source_coord)
        ix = int(round(float(xp))); iy = int(round(float(yp)))
        if 0 <= ix < cube.shape[2] and 0 <= iy < cube.shape[1]:
            center_spec = on_data[:, iy, ix]
            if np.isfinite(center_spec).any():
                peak_center = float(np.nanmax(center_spec))
                snr_center = peak_center / noise if noise > 0 else np.nan

    return {"peak_K": peak, "noise_K": noise,
            "snr": snr_global, "snr_global": snr_global,
            "snr_center": snr_center, "peak_K_center": peak_center,
            "peak_v_kms": v_at_peak}


def make_circle_region(coord, radius_arcsec):
    return regions.CircleSkyRegion(
        coord, radius=radius_arcsec * u.arcsec)


def run_one_source(name, vcen=None, distance=None,
                   region_radius=2.0,
                   families=("NaCl", "KCl", "H2O"),
                   v_upper_max=3,
                   dv_stack=25,
                   make_pv=True,
                   pv_pa=0 * u.deg,
                   pv_length=3 * u.arcsec,
                   include_dihca=True,
                   skip_convolve=False):
    """Run the full salt-search pipeline for *name*.

    Parameters override the source-table defaults so a per-source script can
    say ``vcen=12.5*u.km/u.s`` etc.
    """
    src = sources.get_source(name)
    coord = SkyCoord(src["RAJ2000"] * u.deg, src["DEJ2000"] * u.deg,
                     frame="icrs")
    if vcen is None:
        vcen = sources.get_vlsr(name)
    if distance is None:
        distance = sources.get_distance(name)
    vcen = u.Quantity(vcen, u.km / u.s)
    distance = u.Quantity(distance, u.kpc) if not isinstance(distance, u.Quantity) else distance

    log.info("=== %s ===  ra=%.5f dec=%.5f  v=%s  d=%s",
             name, coord.ra.deg, coord.dec.deg, vcen, distance)

    if np.isnan(vcen.value):
        log.warning("%s: v_LSR unknown — skipping", name)
        skip = paths.source_dir(name) / "report.json"
        with open(skip, "w") as fh:
            json.dump({"source": name, "status": "skipped_no_vlsr"}, fh, indent=2)
        return {"source": name, "status": "skipped_no_vlsr"}

    out = paths.source_dir(name)
    region = make_circle_region(coord, region_radius)
    cube_refs = cube_mod.find_cubes_for_source(
        name, coord.ra.deg, coord.dec.deg, include_dihca=include_dihca)
    log.info("found %d cubes", len(cube_refs))
    if not cube_refs:
        raise FileNotFoundError(
            f"{name}: no cubes under {paths.UVDATA} or {paths.DIHCA_GROUPED}")

    report = {
        "source": name,
        "ra_deg": coord.ra.deg, "dec_deg": coord.dec.deg,
        "vcen_kms": vcen.value, "distance_kpc": distance.value,
        "n_cubes": len(cube_refs),
        "cube_paths": [c.path for c in cube_refs],
        "stacks": {}, "moments": {}, "pv": {},
    }

    for family in families:
        pool = {"NaCl": line_mod.NACL_LINES,
                "KCl":  line_mod.KCL_LINES,
                "H2O":  line_mod.WATER_LINES}[family]
        rests = list(pool.values())
        relevant = []
        for c in cube_refs:
            if any(c.f_lo_GHz <= u.Quantity(r).to(u.GHz).value <= c.f_hi_GHz
                   for r in rests):
                relevant.append(c)
        if not relevant:
            log.info("  %s: no cubes cover any line — skipping", family)
            continue

        log.info("  %s: %d relevant cubes", family, len(relevant))
        try:
            cubes = stacking.open_cubes(relevant, region=region, contsub=True)
            (stack, _, beam), used = stacking.stack_salt_family(
                cubes, family, vcen, dv=dv_stack * u.km / u.s,
                v_upper_max=v_upper_max,
                skip_convolve=skip_convolve,
            )
        except Exception as e:
            log.exception("  %s stacking failed: %s", family, e)
            continue

        stack_path = paths.source_dir(name, "stacks") / f"{family}_stack.fits"
        stacking.write_stack(stack, stack_path)

        # Detection significance: snr_center is the on-axis spectrum at the
        # source pixel; snr_global is the brightest pixel anywhere in the
        # on-velocity slab (kept for diagnostics).
        det = _detection_stats(stack, vcen, source_coord=coord)
        det_flag = (np.isfinite(det["snr_center"]) and det["snr_center"] >= 5)
        report["stacks"][family] = {
            "path": str(stack_path),
            "lines_used": [n for n, _ in used],
            "common_beam": str(beam),
            "peak_K": det["peak_K"],
            "peak_K_center": det["peak_K_center"],
            "noise_K": det["noise_K"],
            "snr": det["snr_global"],
            "snr_global": det["snr_global"],
            "snr_center": det["snr_center"],
            "peak_v_kms": det["peak_v_kms"],
            "detected": bool(det_flag),
        }

        # quicklook moments
        ql_path = paths.source_dir(name, "moments") / f"{family}_quicklook.png"
        try:
            mom_mod.quicklook_grid(stack, ql_path, vcen=vcen,
                                   title=f"{name}: {family} stack")
            report["moments"][family + "_quicklook"] = str(ql_path)
        except Exception as e:
            log.warning("  quicklook failed: %s", e)

        # RGB
        try:
            rgb_path = paths.source_dir(name, "moments") / f"{family}_rgb.png"
            m0b, m0g, m0r = mom_mod.moment_maps(stack, vcen)
            mom_mod.rgb_figure(m0b, m0g, m0r, rgb_path,
                               title=f"{name}: {family} blue/green/red")
            report["moments"][family + "_rgb"] = str(rgb_path)
        except Exception as e:
            log.warning("  RGB failed: %s", e)

        # average spectrum
        try:
            spec_path = paths.source_dir(name, "spectra") / f"{family}_avg_spectrum.png"
            mom_mod.average_spectrum(stack, region=None, outpath=spec_path)
            report["spectra_" + family] = str(spec_path)
        except Exception as e:
            log.warning("  spectrum failed: %s", e)

        # PV
        if make_pv:
            try:
                path = pv_mod.make_path(coord, length=pv_length, pa=pv_pa)
                pv_path = paths.source_dir(name, "pv") / f"{family}_pv.png"
                pv_mod.plot_pv(stack, path, vcen, pv_path,
                               distance=distance,
                               title=f"{name}: {family} PV")
                report["pv"][family] = str(pv_path)
            except Exception as e:
                log.warning("  PV failed: %s", e)

    report_path = paths.source_dir(name) / "report.json"
    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=2, default=str)
    log.info("wrote %s", report_path)
    return report


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("--vcen-kms", type=float, default=None)
    ap.add_argument("--distance-kpc", type=float, default=None)
    ap.add_argument("--region-radius", type=float, default=2.0)
    ap.add_argument("--no-pv", action="store_true")
    ap.add_argument("--no-dihca", action="store_true")
    args = ap.parse_args()
    kw = {}
    if args.vcen_kms is not None:
        kw["vcen"] = args.vcen_kms * u.km / u.s
    if args.distance_kpc is not None:
        kw["distance"] = args.distance_kpc * u.kpc
    run_one_source(args.name, region_radius=args.region_radius,
                   make_pv=not args.no_pv,
                   include_dihca=not args.no_dihca, **kw)
