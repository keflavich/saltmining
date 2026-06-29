"""Per-(target, ALMA program) upper-limit table for the demography_2026 paper.

For each (target, proposal) under analysis_products/, find the brightest mm
continuum source via continuum_sources.csv (max peak_Jybeam, presumed highest
mass), then for each line whose spectrum was recorded under
source_{id:02d}/*.spec.npz, compute two upper limits:

  * UL3_native_mK : 3 sigma peak-T UL at native channel width
  * UL3_10kms_mK  : 3 sigma peak-T UL after smoothing to 10 km/s

For 10 km/s smoothing of a native spectrum at dv_chan,
    sigma_smooth = sigma_native / sqrt(N), with N = max(1, 10/dv_chan)
so a Hanning-smoothed channel UL is 3*sigma_native*sqrt(dv_chan/10).
Lines flagged as detected in line_measurements.csv (snr>=5) are EXCLUDED.

Output:
  data/upper_limits_per_program.csv  -- all rows
  /orange/adamginsburg/salt/demography_2026/upper_limits.tex -- longtable
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
OUT_CSV = ROOT / "data/upper_limits_per_program.csv"
OUT_TEX = Path("/orange/adamginsburg/salt/demography_2026/upper_limits.tex")

KEEP_PREFIXES = ("NaCl_", "KCl_", "H2O_")


def brightest_source_id(cont_csv: Path) -> int | None:
    if not cont_csv.exists() or cont_csv.stat().st_size == 0:
        return None
    try:
        df = pd.read_csv(cont_csv)
    except pd.errors.EmptyDataError:
        return None
    if df.empty or "peak_Jybeam" not in df.columns:
        return None
    return int(df.loc[df["peak_Jybeam"].idxmax(), "id"])


def measure_spec(npz_path: Path) -> dict | None:
    d = np.load(npz_path)
    if not {"vaxis", "spec", "sigma"} <= set(d.files):
        return None
    v = np.asarray(d["vaxis"], dtype=float)
    s = np.asarray(d["spec"], dtype=float)
    sigma_native = float(d["sigma"])
    if v.size < 2 or not np.isfinite(sigma_native):
        return None
    dv = float(np.median(np.abs(np.diff(v))))
    if not np.isfinite(dv) or dv <= 0:
        return None
    N = max(1.0, 10.0 / dv)
    sigma_10 = sigma_native / np.sqrt(N)
    return dict(
        dv_kms=dv,
        sigma_native_K=sigma_native,
        UL3_native_mK=3000.0 * sigma_native,
        UL3_10kms_mK=3000.0 * sigma_10,
        peak_K=float(np.nanmax(s)) if s.size else np.nan,
    )


def collect():
    rows = []
    for proposal_dir in sorted(ANALYSIS.glob("*/2*")):
        target = proposal_dir.parent.name
        proposal = proposal_dir.name
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
            if not line.startswith(KEEP_PREFIXES):
                continue
            if line in detected:
                continue
            m = measure_spec(npz)
            if m is None:
                continue
            rows.append({
                "target": target, "proposal": proposal,
                "brightest_source_id": bid,
                "line": line,
                **m,
            })
    return pd.DataFrame(rows)


def fmt_line(line: str) -> str:
    return line.replace("_", r"\_")


def write_tex(df: pd.DataFrame):
    df = df.sort_values(["target", "proposal", "line"])
    lines = []
    lines.append(r"\begin{longrotatetable}")
    lines.append(r"\begin{deluxetable*}{llllrcc}")
    lines.append(r"\tablecaption{Upper limits on salt (NaCl, KCl) and water "
                 r"emission toward the brightest mm continuum source in each "
                 r"target field, per ALMA program. Detections ($\geq 5\sigma$) "
                 r"are excluded.\label{tab:uls}}")
    lines.append(r"\tablehead{")
    lines.append(r"\colhead{Source} & \colhead{ALMA program} & "
                 r"\colhead{Line} & \colhead{Src ID} & "
                 r"\colhead{$\Delta v_\mathrm{chan}$} & "
                 r"\colhead{$3\sigma$ (native)} & "
                 r"\colhead{$3\sigma$ ($10$\,\kms)} \\")
    lines.append(r"& & & & (\kms) & (mK) & (mK) }")
    lines.append(r"\startdata")
    for _, r in df.iterrows():
        target_esc = r['target'].replace('_', r'\_')
        lines.append(
            f"{target_esc} & {r['proposal']} & "
            f"{fmt_line(r['line'])} & {int(r['brightest_source_id'])} & "
            f"{r['dv_kms']:.2f} & {r['UL3_native_mK']:.1f} & "
            f"{r['UL3_10kms_mK']:.1f} \\\\"
        )
    lines.append(r"\enddata")
    lines.append(r"\tablecomments{The brightest mm continuum source within "
                 r"each (target, ALMA program) is taken as a proxy for the "
                 r"most luminous/massive young stellar object in the field. "
                 r"$3\sigma$ peak-brightness upper limits are reported at the "
                 r"native channel width $\Delta v_\mathrm{chan}$ and after "
                 r"smoothing to a $10$\,\kms\ channel. Smoothing reduces the "
                 r"noise by $\sqrt{\Delta v_\mathrm{chan} / 10\,\mathrm{km\,s^{-1}}}$.}")
    lines.append(r"\end{deluxetable*}")
    lines.append(r"\end{longrotatetable}")
    OUT_TEX.write_text("\n".join(lines) + "\n")


def main():
    df = collect()
    print(f"collected {len(df)} (target, proposal, line) UL rows")
    if df.empty:
        print("nothing to write")
        return
    df.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV}")
    write_tex(df)
    print(f"wrote {OUT_TEX}")
    print()
    print(df.groupby(["target", "proposal"]).size().rename("n_lines"))


if __name__ == "__main__":
    main()
