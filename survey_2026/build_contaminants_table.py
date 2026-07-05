"""LaTeX table of NaCl/KCl/H2O lines vs candidate contaminants.

Reads data/detected_lines_with_contaminants.csv (produced by
find_contaminants_from_xclass.py) and emits a deluxetable to
/orange/adamginsburg/salt/demography_2026/contaminants.tex listing,
for each detected NaCl/KCl/H2O/SiS line, the top contaminant species
from the XCLASS LTE model along with the rest frequency offset between
the contaminant and the observed peak.
"""
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
CSV = ROOT / "data/detected_lines_with_contaminants.csv"
OUT = Path("/orange/adamginsburg/salt/demography_2026/contaminants.tex")

LINE_PATTERN = re.compile(r"^(NaCl|KCl|H2O|SiS)", re.IGNORECASE)


def latex_mol(s):
    s = str(s).strip().rstrip(";")
    # XCLASS labels: CH3OCHO;v=0;ETA1;Sigma1 etc.
    base = s.split(";")[0]
    sub = {
        "CH3OH": r"CH$_3$OH", "CH3CN": r"CH$_3$CN",
        "CH3OCHO": r"CH$_3$OCHO", "CH3OCH3": r"CH$_3$OCH$_3$",
        "C2H5CN": r"C$_2$H$_5$CN", "HC3N": r"HC$_3$N",
        "NH2CHO": r"NH$_2$CHO", "H2CO": r"H$_2$CO",
        "H2O": r"H$_2$O", "SO2": r"SO$_2$",
        "13CO": r"$^{13}$CO", "C18O": r"C$^{18}$O",
        "NaCl": "NaCl", "KCl": "KCl", "SiO": "SiO", "SiS": "SiS",
    }
    rest = s[len(base):].lstrip(";").replace(";", ", ").replace("_", r"\_")
    rest = rest.replace("v=", r"$v$=")
    pretty = sub.get(base, base.replace("_", r"\_"))
    return rf"{pretty}\,{rest}" if rest else pretty


def latex_line(s):
    """Pretty-print the detected line label."""
    s = str(s)
    s = re.sub(r"alpha", r"$\\alpha$", s)
    s = re.sub(r"beta", r"$\\beta$", s)
    s = re.sub(r"gamma", r"$\\gamma$", s)
    s = re.sub(r"H2O", r"H$_2$O", s)
    s = re.sub(r"_v(\d)", r" $v$=\1", s)
    s = re.sub(r"_J(\d+)-(\d+)", r" $J$=\1--\2", s)
    s = re.sub(r"^(NaCl|KCl|SiS|SiO) (.*)", r"\1 \2", s)
    s = s.replace("_", r"\_")
    return s


def main():
    if not CSV.exists():
        print(f"no input file {CSV}")
        return
    df = pd.read_csv(CSV)
    # Filter to NaCl/KCl/H2O/SiS detected lines; one row per (target, line,
    # candidate) with non-zero offset.
    df = df[df["detected_line"].astype(str).apply(lambda s: bool(LINE_PATTERN.match(s)))]
    if df.empty:
        print("no NaCl/KCl/H2O candidate-contaminant matches")
        return
    # Compute |offset|, sort within each (target, line) by smallest |offset|
    # then by largest predicted XCLASS T_K. Take top-3 contaminants per line.
    df = df.assign(abs_off=df["freq_offset_MHz"].abs())
    df = df.sort_values(["target", "detected_line", "abs_off"])
    rows = []
    for (tgt, line), grp in df.groupby(["target", "detected_line"]):
        for _, r in grp.head(3).iterrows():
            rows.append(r)
    if not rows:
        return
    table = pd.DataFrame(rows)
    out = []
    out.append(r"\startlongtable")
    out.append(r"\begin{deluxetable}{llcccc}")
    out.append(r"\tabletypesize{\scriptsize}")
    out.append(r"\tablecaption{Candidate contaminants for NaCl, KCl, "
               r"H$_2$O, and SiS detections. For every detected target "
               r"line in these species at the brightest mm continuum "
               r"source we list up to three nearest-in-frequency "
               r"transitions returned by the XCLASS LTE model "
               r"(\texttt{transition\_energies.dat}); the offset is "
               r"between the candidate's source-frame Doppler-shifted "
               r"frequency and the observed line peak. Predicted "
               r"$T_\mathrm{B}$ is the candidate's modeled XCLASS "
               r"brightness temperature, included as a guide to which "
               r"contaminants the LTE model considers strongest in "
               r"the band.\label{tab:contaminants}}")
    out.append(r"\tablehead{")
    out.append(r"\colhead{Source} & \colhead{Detected line} & "
               r"\colhead{Rest freq.} & \colhead{Contaminant} & "
               r"\colhead{$\Delta\nu$} & \colhead{$T_\mathrm{B}^{XCLASS}$} \\")
    out.append(r" & & (GHz) & & (MHz) & (K) }")
    out.append(r"\startdata")
    cur_target = None
    for _, r in table.iterrows():
        tgt = str(r["target"])
        line = latex_line(r["detected_line"])
        rest_g = float(r["detected_rest_GHz"])
        cand = latex_mol(r["candidate_mol"])
        off = float(r["freq_offset_MHz"])
        tb = float(r["candidate_T_K"])
        src_cell = tgt.replace("_", r"\_") if tgt != cur_target else ""
        cur_target = tgt
        out.append(f"{src_cell} & {line} & {rest_g:.4f} & {cand} & "
                   f"{off:+.2f} & {tb:.1f} \\\\")
    out.append(r"\enddata")
    out.append(r"\tablecomments{Source positions and XCLASS LTE model "
               r"parameters as listed in \texttt{lineid/<target>/}. "
               r"The contaminant catalog spans only the species and "
               r"vibrational states enabled in the model "
               r"(\texttt{<target>\_band6.molfit}); blends with "
               r"transitions outside that list are not captured. "
               r"Self-matches within $\pm5$\,MHz of the detected line "
               r"have been suppressed.}")
    out.append(r"\end{deluxetable}")
    OUT.write_text("\n".join(out) + "\n")
    print(f"wrote {OUT} ({len(table)} rows, "
          f"{table['target'].nunique()} sources)")


if __name__ == "__main__":
    main()
