"""Per-source line inventory + contaminant flags + probable-detection score.

For every (target, proposal, brightest_src) in analysis_products that has
a populated line_measurements.csv, emit:

  data/per_source_linelist.csv      machine-readable per-line table
  /orange/adamginsburg/salt/demography_2026/per_source_linelist.tex
                                       human-readable LaTeX table grouped
                                       by source

Species included (target list): NaCl, KCl, SiS, SiO, H2O, H2S, PN, RRL.
Context species (logged but not asserted as detections): SO, COMs.

A "probable detection" requires:
  - snr >= 5
  - peak_v within ±5 km/s of source-systemic (when known)
  - no XCLASS contaminant within ±5 MHz whose predicted T_B exceeds
    0.3 * observed peak T
  - if vibrationally-excited and v=0 is in band: v=0 must also be detected
    (else flagged as vexc-only suspect, dropped from probable list)
"""
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
CONT_CSV = ROOT / "data/detected_lines_with_contaminants.csv"
VLSR_LIT = ROOT / "data/vlsr_from_literature.json"
VLSR_DATA = ROOT / "data/vlsr_from_data.json"
OUT_CSV = ROOT / "data/per_source_linelist.csv"
OUT_TEX = Path("/orange/adamginsburg/salt/demography_2026/per_source_linelist.tex")

SPECIES_GROUPS = {
    "NaCl": re.compile(r"^Na(37)?Cl(_|v=)"),
    "KCl":  re.compile(r"^(41)?K(37)?Cl(_|v=)"),
    "SiS":  re.compile(r"^SiS"),
    "SiO":  re.compile(r"^Si(29|30)?O(_|v)"),
    "H2O":  re.compile(r"^H2O"),
    "H2S":  re.compile(r"^H2S"),
    "PN":   re.compile(r"^PN(_|v)"),
    "RRL":  re.compile(r"^H\d+(alpha|beta|gamma|delta)$"),
    "SO":   re.compile(r"^SO[_0-9]"),
    "COM":  re.compile(r"^(CH3OH|CH3CN|CH3OCHO|CH3OCH3|C2H5CN|HC3N|NH2CHO|HCOOH|HNCO|H2CO)"),
}
SCORE_SPECIES = ("NaCl", "KCl", "SiS", "H2O", "H2S", "PN", "RRL")

C_KMS = 299792.458


def species_of(line):
    s = str(line)
    for sp, rx in SPECIES_GROUPS.items():
        if rx.match(s):
            return sp
    return None


def vexc(line):
    m = re.search(r"_v(\d+)_|v=(\d+)_", str(line))
    if m:
        return int(m.group(1) or m.group(2))
    return -1


def brightest_id(prop_dir):
    cont = prop_dir / "continuum_sources.csv"
    if not cont.exists() or cont.stat().st_size == 0:
        return None
    try:
        df = pd.read_csv(cont)
    except pd.errors.EmptyDataError:
        return None
    if df.empty or "peak_Jybeam" not in df.columns:
        return None
    return int(df.loc[df["peak_Jybeam"].idxmax(), "id"])


def load_vlsr_lookup():
    import json
    out = {}
    for p in (VLSR_LIT, VLSR_DATA):
        if p.exists():
            try:
                d = json.loads(p.read_text())
                for k, v in d.items():
                    val = v.get("v_LSR_kms")
                    if val is not None:
                        out[k] = float(val)
            except (json.JSONDecodeError, AttributeError):
                pass
    return out


def load_contaminants():
    if not CONT_CSV.exists():
        return pd.DataFrame()
    return pd.read_csv(CONT_CSV)


def contaminant_for(cont_df, target, line, max_offset_MHz=5.0):
    """Return (top_candidate, best_offset_MHz, top_TB) or (None, None, None)
    for a (target, line). Picks the candidate with the smallest |Δν|."""
    if cont_df.empty:
        return None, None, None
    sel = cont_df[(cont_df["target"] == target)
                   & (cont_df["detected_line"] == line)]
    if sel.empty:
        return None, None, None
    sel = sel.assign(abs_off=sel["freq_offset_MHz"].abs())
    sel = sel.sort_values("abs_off")
    best = sel.iloc[0]
    if abs(float(best["freq_offset_MHz"])) > max_offset_MHz:
        return None, None, None
    return (str(best["candidate_mol"]).strip(),
            float(best["freq_offset_MHz"]),
            float(best["candidate_T_K"]))


def main():
    cont_df = load_contaminants()
    vlsr_lookup = load_vlsr_lookup()
    rows = []
    for target_dir in sorted(ANALYSIS.iterdir()):
        if not target_dir.is_dir():
            continue
        target = target_dir.name
        for prop_dir in sorted(target_dir.glob("2*")):
            meas = prop_dir / "line_measurements.csv"
            if not meas.exists():
                continue
            try:
                df = pd.read_csv(meas)
            except pd.errors.EmptyDataError:
                continue
            if df.empty or "snr" not in df.columns:
                continue
            bid = brightest_id(prop_dir)
            if bid is None:
                continue
            sub = df[(df["source"] == bid)
                      & (df["snr"] >= 5.0)].copy()
            if sub.empty:
                continue
            v_sys = vlsr_lookup.get(target)
            # Pre-compute v=0 presence per species for this proposal
            v0_present = {}
            for _, r in sub.iterrows():
                sp = species_of(r["line"])
                if sp and vexc(r["line"]) == 0:
                    v0_present[sp] = True
            for _, r in sub.iterrows():
                line = str(r["line"])
                sp = species_of(line)
                if sp is None:
                    continue
                snr = float(r["snr"])
                rest_GHz = float(r["rest_GHz"])
                peak_v = (float(r["peak_v"])
                          if pd.notna(r.get("peak_v")) else np.nan)
                peak_T = float(r.get("peak_Kkms_or_unit", np.nan))
                vex = vexc(line)
                # Contaminant flag
                top_cand, top_off, top_TB = contaminant_for(
                    cont_df, target, line, max_offset_MHz=5.0)
                contam_severity = "low"
                if top_cand is not None:
                    if not np.isfinite(peak_T) or peak_T <= 0:
                        contam_severity = "med"
                    elif top_TB >= 0.3 * abs(peak_T):
                        contam_severity = "high"
                    else:
                        contam_severity = "med"
                # Velocity offset flag
                if pd.notna(peak_v) and v_sys is not None:
                    dv = peak_v - v_sys
                    v_off_flag = abs(dv) > 5.0
                else:
                    dv = np.nan
                    v_off_flag = False
                # Vexc-only flag (v>=2 detected, v=0 in band but undetected)
                vexc_only = False
                if vex >= 2:
                    if not v0_present.get(sp, False):
                        vexc_only = True
                # Probable detection
                probable = (snr >= 5.0
                             and sp in SCORE_SPECIES
                             and contam_severity != "high"
                             and not v_off_flag
                             and not vexc_only)
                rows.append(dict(
                    target=target, proposal=prop_dir.name, src=bid,
                    species=sp, line=line, rest_GHz=rest_GHz, snr=snr,
                    peak_T=peak_T, peak_v=peak_v, dv=dv, v_state=vex,
                    contaminant=top_cand if top_cand else "",
                    contam_offset_MHz=top_off if top_off is not None else "",
                    contam_TB=top_TB if top_TB is not None else "",
                    contam_severity=contam_severity,
                    vexc_only=vexc_only,
                    v_off_flag=v_off_flag,
                    probable=probable,
                ))
    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV} ({len(out)} rows, "
          f"{out['target'].nunique() if not out.empty else 0} sources)")

    # Build LaTeX table grouped by source
    lines = []
    lines.append(r"\startlongtable")
    lines.append(r"\begin{deluxetable}{lllcccll}")
    lines.append(r"\tabletypesize{\scriptsize}")
    lines.append(r"\tablecaption{Complete per-source line inventory for the "
                 r"NaCl, KCl, SiS, SiO, H$_2$O, H$_2$S, PN, and RRL search "
                 r"set. Each row is a $\geq 5\sigma$ line at the brightest "
                 r"mm continuum source of the listed proposal. `probable' "
                 r"detections satisfy: low XCLASS-contaminant severity "
                 r"($T_B^{\rm cand}<0.3\,T_B^{\rm obs}$ or no candidate "
                 r"within $\pm5$\,MHz), velocity within $\pm5$\,\kms\ of "
                 r"the source $v_\mathrm{LSR}$ when known, and v=0 "
                 r"corroboration if a v$\geq$2 line is the detection."
                 r"\label{tab:per_source_linelist}}")
    lines.append(r"\tablehead{")
    lines.append(r"\colhead{Source} & \colhead{Species} & "
                 r"\colhead{Line} & \colhead{$\nu_0$ (GHz)} & "
                 r"\colhead{SNR} & \colhead{$T_B$ (K)} & "
                 r"\colhead{Contam.} & \colhead{Flag}}")
    lines.append(r"\startdata")
    try:
        import build_nacl_rrl_table as v1
    except ImportError:
        v1 = None
    last_target = None
    for _, r in (out.sort_values(["target", "species", "rest_GHz"])
                  .iterrows() if not out.empty else iter([])):
        raw_target = str(r["target"])
        target = v1._display_name(raw_target) if v1 else raw_target
        first_row = (target != last_target)
        last_target = target
        line_tex = (str(r["line"]).replace("_", r"\_"))
        contam = str(r.get("contaminant", "") or "")
        contam_tex = contam.replace("_", r"\_") if contam else r"\nodata"
        flags = []
        if r["probable"]:
            flags.append(r"\textbf{prob}")
        else:
            if r["contam_severity"] == "high":
                flags.append("blend")
            if r["vexc_only"]:
                flags.append("vexc-only")
            if r["v_off_flag"]:
                flags.append(rf"$\Delta v$={r['dv']:.1f}")
        flag_str = ", ".join(flags) if flags else ""
        src_cell = target.replace("_", r"\_") if first_row else ""
        T_str = (f"{r['peak_T']:.3f}"
                  if isinstance(r['peak_T'], (int, float))
                     and np.isfinite(r['peak_T']) else r"\nodata")
        lines.append(
            f"{src_cell} & {r['species']} & {line_tex} & "
            f"{r['rest_GHz']:.4f} & {r['snr']:.1f} & "
            f"{T_str} & {contam_tex} & {flag_str} \\\\"
        )
    lines.append(r"\enddata")
    lines.append(r"\end{deluxetable}")
    OUT_TEX.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT_TEX}")

    # Summary count
    if not out.empty:
        print()
        prob = out[out["probable"]]
        for sp in SCORE_SPECIES:
            n_total = (out["species"] == sp).sum()
            n_prob = (prob["species"] == sp).sum() if not prob.empty else 0
            print(f"  {sp:5s}: {n_total:4d} >=5sig rows, "
                  f"{n_prob:4d} probable")


if __name__ == "__main__":
    main()
