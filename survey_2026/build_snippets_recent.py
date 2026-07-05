"""Build CARTA snippets for recently-downloaded regions.

For each target, find local cont + cube files across all downloaded proposals
and write a snippet to ~/.carta/config/snippets/ under category saltsearch2026.
"""
import json
import re
from pathlib import Path

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
UVDIR = ROOT / "uvdata"
SNIP = Path.home() / ".carta/config/snippets"
SNIP.mkdir(parents=True, exist_ok=True)

CATEGORY = "saltsearch2026"

TARGETS = None  # set None to iterate every L4_d2 target with downloaded data


def normalize(s):
    return re.sub(r"[\s_\-+]", "", str(s)).lower()


def name_variants(s):
    base = str(s)
    out = {normalize(base),
           normalize(base.replace("+", "p").replace("-", "m")),
           normalize(base.replace("+", "").replace("-", ""))}
    return {v for v in out if v}


def _brightest_source_coord(target):
    """Return (ra_deg, dec_deg) of the brightest mm continuum source across
    all analyzed proposals, or None."""
    import pandas as pd
    tgt_dir = ROOT / "analysis_products" / target
    if not tgt_dir.is_dir():
        return None
    best = None
    for prop_dir in sorted(tgt_dir.glob("2*")):
        cont = prop_dir / "continuum_sources.csv"
        if not cont.exists():
            continue
        try:
            df = pd.read_csv(cont)
        except pd.errors.EmptyDataError:
            continue
        if df.empty or "peak_Jybeam" not in df.columns:
            continue
        idx = int(df["peak_Jybeam"].idxmax())
        peak = float(df.loc[idx, "peak_Jybeam"])
        ra = float(df.loc[idx, "ra_deg"])
        dec = float(df.loc[idx, "dec_deg"])
        if best is None or peak > best[0]:
            best = (peak, ra, dec)
    if best is None:
        return None
    return (best[1], best[2])


def _cube_in_fov(cube_path, ra, dec):
    """True if the (ra, dec) world coord falls inside cube_path's FOV."""
    from astropy.io import fits
    from astropy.wcs import WCS
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter("ignore")
        h = fits.getheader(cube_path)
        try:
            wcs = WCS(h).celestial
            xp, yp = wcs.world_to_pixel_values(ra, dec)
        except (ValueError, KeyError):
            return False
        nx = int(h.get("NAXIS1", 0))
        ny = int(h.get("NAXIS2", 0))
        return (0 <= float(xp) < nx) and (0 <= float(yp) < ny)


_CAL_PATTERNS = (
    re.compile(r"J\d{4}[-+]\d{4}"),         # quasar calibrators
    re.compile(r"_ph\.", re.IGNORECASE),    # phase cal
    re.compile(r"_ampph\.", re.IGNORECASE),
    re.compile(r"_bp\.", re.IGNORECASE),    # bandpass cal
    re.compile(r"_check\.", re.IGNORECASE),
    re.compile(r"_pol\.", re.IGNORECASE),
)


def _is_calibrator(filename: str) -> bool:
    return any(p.search(filename) for p in _CAL_PATTERNS)


def find_files(target):
    """Return (best_cont, sorted_cubes) across all proposals for this target.
    Filters out cubes whose FOV does NOT contain the brightest mm continuum
    source for this target (handles multi-pointing projects where the same
    directory holds data for unrelated fields), and skips ALMA calibrator
    files (quasar phase/bandpass/check cals)."""
    tgt_keys = name_variants(target)
    coord = _brightest_source_coord(target)
    cont_paths, cube_paths = [], []
    for prop_dir in sorted(UVDIR.iterdir()):
        if not prop_dir.is_dir():
            continue
        for sub in prop_dir.iterdir():
            if not sub.is_dir():
                continue
            if normalize(sub.name) not in tgt_keys:
                continue
            for f in sub.iterdir():
                if not f.name.endswith(".I.pbcor.fits"):
                    continue
                if _is_calibrator(f.name):
                    continue
                if ".cont." in f.name:
                    cont_paths.append(f)
                elif ".cube." in f.name:
                    cube_paths.append(f)
    # FOV filter
    if coord is not None:
        ra, dec = coord
        cont_paths = [p for p in cont_paths if _cube_in_fov(p, ra, dec)]
        cube_paths = [p for p in cube_paths if _cube_in_fov(p, ra, dec)]
    cont_best = None
    if cont_paths:
        multi = [p for p in cont_paths if re.search(r"spw\d+(_\d+){2,}", p.name)]
        pool = multi if multi else cont_paths
        cont_best = sorted(pool, key=lambda p: p.stat().st_size)[0]
    return cont_best, sorted(cube_paths, key=lambda p: p.name)


def find_field_stack_fits(target):
    """Field-level stack cube FITS files (NaCl/KCl/H2O) from analysis_products/
    <target>/stacks/. These are 3D cubes that CARTA can display."""
    sdir = ROOT / "analysis_products" / target / "stacks"
    if not sdir.is_dir():
        return []
    return sorted(sdir.glob("*_stack.fits"))


CARTA_ROOT_PREFIX = "/orange/adamginsburg"


def carta_path(p):
    s = str(p)
    return s[len(CARTA_ROOT_PREFIX):] if s.startswith(CARTA_ROOT_PREFIX) else s


def find_region_files(target):
    """Return sorted list of *.reg files under analysis_products/<target>/<prop>/."""
    base = ROOT / "analysis_products" / target
    if not base.is_dir():
        return []
    regs = []
    for prop_dir in sorted(base.glob("2*")):
        for r in prop_dir.glob("*.reg"):
            regs.append(r)
    return regs


def find_detection_mom0s(target, max_per_source=3, snr_threshold=5.0):
    """For each (proposal, source) with >=1 line at snr>=5, return up to
    max_per_source mom0 FITS paths (highest-SNR lines first). New-naming
    convention: source_<NN>_<line>_mom0.fits inside source_<NN>/."""
    import pandas as pd
    base = ROOT / "analysis_products" / target
    if not base.is_dir():
        return []
    out = []
    for prop_dir in sorted(base.glob("2*")):
        meas = prop_dir / "line_measurements.csv"
        if not meas.exists():
            continue
        try:
            df = pd.read_csv(meas)
        except pd.errors.EmptyDataError:
            continue
        if df.empty or "snr" not in df.columns:
            continue
        hits = df[df["snr"] >= snr_threshold]
        if hits.empty:
            continue
        for sid, g in hits.groupby("source"):
            sid_i = int(sid)
            top = g.nlargest(max_per_source, "snr")
            sdir = prop_dir / f"source_{sid_i:02d}"
            for _, r in top.iterrows():
                p_new = sdir / f"source_{sid_i:02d}_{r['line']}_mom0.fits"
                p_old = sdir / f"{r['line']}_mom0.fits"
                if p_new.exists():
                    out.append(p_new)
                elif p_old.exists():
                    out.append(p_old)
    return out


def build_snippet_code(cont, cubes, regs, mom0s=None, stacks=None):
    lines = []
    rest = list(cubes)
    if cont is not None:
        lines.append(f'await app.openFile("{carta_path(cont)}")')
    elif rest:
        first = rest.pop(0)
        lines.append(f'await app.openFile("{carta_path(first)}")')
    for c in rest:
        lines.append(f'await app.appendFile("{carta_path(c)}")')
    if mom0s:
        for m in mom0s:
            lines.append(f'await app.appendFile("{carta_path(m)}")')
    if stacks:
        for s in stacks:
            lines.append(f'await app.appendFile("{carta_path(s)}")')
    for r in regs:
        rdir = carta_path(r.parent) + "/"
        rname = r.name
        lines.append(f'await app.importRegion("{rdir}", "{rname}", 2)')
    return "\n".join(lines)


def main():
    if TARGETS is None:
        import pandas as pd
        df = pd.read_csv(ROOT / "data/sources_L4_d2.csv")
        targets = sorted(set(df["name"]))
    else:
        targets = TARGETS
    for tgt in targets:
        cont, cubes = find_files(tgt)
        regs = find_region_files(tgt)
        if cont is None and not cubes:
            print(f"  {tgt}: NO LOCAL DATA — skipping")
            continue
        mom0s = find_detection_mom0s(tgt)
        stacks = find_field_stack_fits(tgt)
        code = build_snippet_code(cont, cubes, regs, mom0s=mom0s,
                                    stacks=stacks)
        snippet = {
            "$schema": "https://cartavis.github.io/schemas/snippet_schema_1.json",
            "categories": [CATEGORY],
            "code": code,
            "frontendVersion": "3.0.0",
            "snippetVersion": 1,
        }
        out = SNIP / f"{tgt}_2026.json"
        out.write_text(json.dumps(snippet, indent=4))
        print(f"  wrote {out.name}: cont={'yes' if cont else 'no'}, "
              f"cubes={len(cubes)}, regs={len(regs)}, mom0s={len(mom0s)}, "
              f"stacks={len(stacks)}")


if __name__ == "__main__":
    main()
