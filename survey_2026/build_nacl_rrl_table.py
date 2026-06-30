"""NaCl + RRL detection/UL table.

For each (target, ALMA program) under analysis_products/, identify the brightest
mm continuum source (proxy for the most massive YSO), then for every NaCl_*
and RRL (H{n}alpha/beta/gamma/delta) line whose spectrum was recorded under
source_{id:02d}/*.spec.npz, report either a detection peak-T or a 3 sigma
upper limit under four conditions:

  T_1kms_native  : 1 km/s channel at the native synthesized beam
  T_10kms_native : 10 km/s channel at the native synthesized beam
  T_1kms_300au   : 1 km/s channel at a 300-AU effective beam (spatially smoothed)
  T_10kms_300au  : 10 km/s channel at a 300-AU effective beam

Scalings used:
  sigma_T(dv) = sigma_T(dv_native) * sqrt(dv_native / dv)        (for dv >= dv_native)
  sigma_T(B_300au) = sigma_T(B_native) * (B_native / B_300)^2    (point-source
                                                                  brightness dilution
                                                                  when smoothing)
  When B_native >= 300 AU at the source distance, the 300-AU column reports
  the native value (cannot resolve finer than native).
  When dv_native >= dv (target finer than native), reports native.

Detections (snr >= 5 in line_measurements.csv) get the peak-T value; rest get
3 sigma ULs. All values in K.

Output:
  data/nacl_rrl_table.csv  (machine-readable)
  /orange/adamginsburg/salt/demography_2026/nacl_rrl_uls.tex  (longtable)
"""
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
SRC_CSV = ROOT / "data/sources_L4_d2.csv"
OUT_CSV = ROOT / "data/nacl_rrl_table.csv"
OUT_TEX = Path("/orange/adamginsburg/salt/demography_2026/nacl_rrl_uls.tex")

RRL_RE = re.compile(r"^H\d+(alpha|beta|gamma|delta)$")


def is_keep_line(line: str) -> bool:
    return line.startswith("NaCl_") or bool(RRL_RE.match(line))


def brightest_source_id(cont_csv: Path):
    if not cont_csv.exists() or cont_csv.stat().st_size == 0:
        return None
    try:
        df = pd.read_csv(cont_csv)
    except pd.errors.EmptyDataError:
        return None
    if df.empty or "peak_Jybeam" not in df.columns:
        return None
    return int(df.loc[df["peak_Jybeam"].idxmax(), "id"])


_BEAM_CACHE: dict[str, float] = {}


def beam_arcsec_for_cube(cube_path: str) -> float | None:
    """Return the synthesized-beam major axis in arcsec for a cube, cached."""
    if cube_path in _BEAM_CACHE:
        return _BEAM_CACHE[cube_path]
    p = Path(cube_path)
    if not p.exists():
        _BEAM_CACHE[cube_path] = None
        return None
    try:
        h = fits.getheader(p, ext=0)
    except (OSError, fits.VerifyError):
        _BEAM_CACHE[cube_path] = None
        return None
    bmaj = h.get("BMAJ")
    if bmaj is None:
        _BEAM_CACHE[cube_path] = None
        return None
    _BEAM_CACHE[cube_path] = float(bmaj) * 3600.0  # deg -> arcsec
    return _BEAM_CACHE[cube_path]


def beam_for_proposal(proposal_dir: Path) -> float | None:
    """Estimate a representative beam (arcsec) from the smallest BMAJ across
    cube files in this proposal directory."""
    # try first cube; the line_pipeline already chose one beam per cube
    cubes = sorted(proposal_dir.glob("../../uvdata/*/*/member.*sci.*cube*pbcor.fits"))
    if not cubes:
        return None
    beams = []
    for c in cubes:
        b = beam_arcsec_for_cube(str(c))
        if b is not None and np.isfinite(b):
            beams.append(b)
    return float(np.min(beams)) if beams else None


def measure_native(npz_path: Path):
    d = np.load(npz_path)
    if not {"vaxis", "spec", "sigma"} <= set(d.files):
        return None
    v = np.asarray(d["vaxis"], dtype=float)
    s = np.asarray(d["spec"], dtype=float)
    sigma = float(d["sigma"])
    if v.size < 2 or not np.isfinite(sigma):
        return None
    dv = float(np.median(np.abs(np.diff(v))))
    if not np.isfinite(dv) or dv <= 0:
        return None
    return dict(dv_native=dv, sigma_native=sigma, peak=float(np.nanmax(s)))


_ALT_BY_TARGET = None


def _display_name(name: str) -> str:
    """Match the naming used by Tables 1/4/5 (COMMON_NAME → IRAS → target)."""
    global _ALT_BY_TARGET
    if _ALT_BY_TARGET is None:
        csv = ROOT / "data/sources_L4_d2.csv"
        if csv.exists():
            try:
                df = pd.read_csv(csv)
                col = df["alma_target_names"] if "alma_target_names" in df.columns else pd.Series([""] * len(df))
                _ALT_BY_TARGET = dict(zip(df["name"], col.fillna("")))
            except pd.errors.EmptyDataError:
                _ALT_BY_TARGET = {}
        else:
            _ALT_BY_TARGET = {}
    try:
        from build_target_table import COMMON_NAME, alma_target_names_to_iras
    except ImportError:
        COMMON_NAME = {}
        alma_target_names_to_iras = lambda _s: None
    base = COMMON_NAME.get(name)
    if not base:
        base = alma_target_names_to_iras(_ALT_BY_TARGET.get(name, "")) or name
    return base


def sigma_at_dv(sigma_native: float, dv_native: float, dv_target: float) -> float:
    if dv_target <= dv_native:
        return sigma_native
    return sigma_native * np.sqrt(dv_native / dv_target)


def sigma_at_300au(sigma: float, beam_native_au: float, target_au: float = 300.0) -> float:
    if not np.isfinite(beam_native_au) or beam_native_au >= target_au:
        return sigma
    return sigma * (beam_native_au / target_au) ** 2


def fmt_value(peak_K, sigma_K, detected):
    """Detection -> peak (mK); UL -> '<' 3sigma (mK)."""
    if detected:
        return f"{peak_K*1000.0:.1f}"
    return f"$<${3.0*sigma_K*1000.0:.1f}"


def collect(distances: dict[str, float]):
    rows = []
    for proposal_dir in sorted(ANALYSIS.glob("*/2*")):
        target = proposal_dir.parent.name
        proposal = proposal_dir.name
        d_kpc = distances.get(target)
        if d_kpc is None or not np.isfinite(d_kpc):
            continue
        cont = proposal_dir / "continuum_sources.csv"
        bid = brightest_source_id(cont)
        if bid is None:
            continue
        src_dir = proposal_dir / f"source_{bid:02d}"
        if not src_dir.is_dir():
            continue
        meas_csv = proposal_dir / "line_measurements.csv"
        detected = set()
        if meas_csv.exists():
            mdf = pd.read_csv(meas_csv)
            if "snr" in mdf.columns:
                hit = mdf[(mdf["source"] == bid) & (mdf["snr"] >= 5.0)]
                detected = set(hit["line"].astype(str))
        for npz in sorted(src_dir.glob("*.spec.npz")):
            line = npz.name.replace(".spec.npz", "")
            if not is_keep_line(line):
                continue
            m = measure_native(npz)
            if m is None:
                continue
            # native beam from continuum source's spectral cube (proxy: any cube)
            beam_arcsec = None
            if meas_csv.exists():
                try:
                    mdf = pd.read_csv(meas_csv)
                    cubes_for_line = mdf[mdf["line"] == line]["cube"].dropna().tolist()
                    if cubes_for_line:
                        cube_full = ROOT / "uvdata" / proposal / target / cubes_for_line[0]
                        beam_arcsec = beam_arcsec_for_cube(str(cube_full))
                except (KeyError, pd.errors.EmptyDataError):
                    pass
            beam_au = (beam_arcsec * d_kpc * 1000.0) if beam_arcsec else np.nan
            sn = m["sigma_native"]
            dvn = m["dv_native"]
            sig_1 = sigma_at_dv(sn, dvn, 1.0)
            sig_10 = sigma_at_dv(sn, dvn, 10.0)
            sig_1_300 = sigma_at_300au(sig_1, beam_au, 300.0)
            sig_10_300 = sigma_at_300au(sig_10, beam_au, 300.0)
            det = line in detected
            rows.append({
                "target": target, "proposal": proposal, "line": line,
                "src_id": bid, "dv_native_kms": dvn,
                "beam_native_arcsec": beam_arcsec,
                "beam_native_au": beam_au,
                "peak_K": m["peak"], "detected": det,
                "sigma_1kms_native_K": sig_1,
                "sigma_10kms_native_K": sig_10,
                "sigma_1kms_300au_K": sig_1_300,
                "sigma_10kms_300au_K": sig_10_300,
            })
    return pd.DataFrame(rows)


def fmt_line(line):
    return line.replace("_", r"\_")


def write_tex(df: pd.DataFrame):
    df = df.sort_values(["target", "proposal", "line"])
    out = []
    out.append(r"\startlongtable")
    out.append(r"\begin{deluxetable*}{lllcccccc}")
    out.append(r"\rotate")
    out.append(r"\tablecaption{NaCl and radio recombination-line (RRL) detections and "
               r"$3\sigma$ upper limits toward the brightest mm continuum source in "
               r"each target field, per ALMA program. Values in mK; \nodata\ indicates "
               r"that the corresponding line was not in the recorded spectrum set. "
               r"Native columns use the cube's own synthesized beam; the $300$\,AU "
               r"columns smooth the cube spatially to a $300$\,AU effective beam at "
               r"the source distance (point-source brightness scaling).\label{tab:naclrrl}}")
    out.append(r"\tablehead{")
    out.append(r"\colhead{Source} & \colhead{Program} & \colhead{Line} & "
               r"\colhead{Src} & \colhead{$\theta_\mathrm{beam}$} & "
               r"\colhead{$1$\,\kms, nat.} & \colhead{$10$\,\kms, nat.} & "
               r"\colhead{$1$\,\kms, $300$\,AU} & \colhead{$10$\,\kms, $300$\,AU} \\")
    out.append(r" & & & & (AU) & (mK) & (mK) & (mK) & (mK) }")
    out.append(r"\startdata")
    for _, r in df.iterrows():
        det = bool(r["detected"])
        v1n = fmt_value(r["peak_K"], r["sigma_1kms_native_K"], det)
        v10n = fmt_value(r["peak_K"], r["sigma_10kms_native_K"], det)
        v1_3 = fmt_value(r["peak_K"], r["sigma_1kms_300au_K"], det)
        v10_3 = fmt_value(r["peak_K"], r["sigma_10kms_300au_K"], det)
        beam_au = r["beam_native_au"]
        beam_str = f"{beam_au:.0f}" if np.isfinite(beam_au) else r"\nodata"
        target = _display_name(str(r["target"])).replace("_", r"\_")
        out.append(
            f"{target} & {r['proposal']} & {fmt_line(r['line'])} & "
            f"{int(r['src_id'])} & {beam_str} & "
            f"{v1n} & {v10n} & {v1_3} & {v10_3} \\\\"
        )
    out.append(r"\enddata")
    out.append(r"\tablecomments{For each row, the brightest mm continuum source in the "
               r"target field is taken as a proxy for the most massive YSO. "
               r"Detections ($\mathrm{S/N} \geq 5$) report the peak brightness temperature; "
               r"non-detections report $3\sigma$ upper limits. Channel-width scaling: "
               r"$\sigma_T(\Delta v) = \sigma_T(\Delta v_\mathrm{nat}) \sqrt{\Delta v_\mathrm{nat}/\Delta v}$. "
               r"Spatial-smoothing scaling for an unresolved disk: "
               r"$\sigma_T(\theta_{300}) = \sigma_T(\theta_\mathrm{nat}) "
               r"(\theta_\mathrm{nat}/\theta_{300})^2$. "
               r"When $\theta_\mathrm{nat} \geq 300$\,AU the native value is reported in the "
               r"$300$\,AU column.}")
    out.append(r"\end{deluxetable*}")
    OUT_TEX.write_text("\n".join(out) + "\n")


def main():
    src = pd.read_csv(SRC_CSV)
    distances = dict(zip(src["name"], src["dist_kpc"]))
    df = collect(distances)
    print(f"collected {len(df)} rows")
    if df.empty:
        print("nothing")
        return
    df.to_csv(OUT_CSV, index=False)
    write_tex(df)
    print(f"wrote {OUT_CSV}")
    print(f"wrote {OUT_TEX}")
    print()
    print(df.groupby(["target", "proposal"]).size().rename("n_lines"))


if __name__ == "__main__":
    main()
