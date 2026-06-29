"""Build the sample-summary table:

  rows = {d <= 2 kpc, d <= 3 kpc}
  cols = N, N_with_ALMA, Analyzed, RRL, NaCl, KCl, H2O, COM, SiS

Detection counts are evaluated AT THE BRIGHTEST mm continuum source of each
target, mirroring build_data_summary. Literature-confirmed detections
(data/literature_detections.csv) take precedence over pipeline results.

Per-species cell: "{det_pct}% ({N_det}/{N_obs})", where N_obs counts targets
in which the line falls within spectral coverage of at least one analyzed
proposal. Sources without spectral coverage of a line are excluded from that
species' denominator.
"""
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
L4_CSV = ROOT / "data/sources_L4_d2.csv"
LIT_CSV = ROOT / "data/literature_detections.csv"
MD3 = Path("/orange/adamginsburg/galactic_plane_surveys/uchii_2026/sample_lt3kpc.md")
OUT_TEX = Path("/orange/adamginsburg/salt/demography_2026/sample_summary.tex")

RRL_RE = re.compile(r"^H\d+(alpha|beta|gamma|delta)$")
SIS_RE = re.compile(r"^SiS", re.IGNORECASE)
COM_PAT = re.compile(r"^(CH3OH|CH3CN|HC3N|CH3OCHO|C2H5CN|CH3OCH3|NH2CHO|HCOOH|HNCO)",
                     re.IGNORECASE)

SPECIES = ("rrl", "nacl", "kcl", "h2o", "com", "sis")


def brightest_id(proposal_dir):
    cont = proposal_dir / "continuum_sources.csv"
    if not cont.exists() or cont.stat().st_size == 0:
        return None
    try:
        df = pd.read_csv(cont)
    except pd.errors.EmptyDataError:
        return None
    if df.empty or "peak_Jybeam" not in df.columns:
        return None
    return int(df.loc[df["peak_Jybeam"].idxmax(), "id"])


def species_of(line, group=None):
    s = str(line)
    g = str(group or "")
    if g == "NaCl":
        return "nacl"
    if g == "KCl":
        return "kcl"
    if g == "H2O":
        return "h2o"
    if g == "RRL" or RRL_RE.match(s):
        return "rrl"
    if SIS_RE.match(s):
        return "sis"
    if COM_PAT.match(s):
        return "com"
    return None


def load_lit():
    """Return dict target -> {species: 'det'|'tent'|'ul'|'na'}."""
    if not LIT_CSV.exists():
        return {}
    df = pd.read_csv(LIT_CSV)
    out = {}
    for _, r in df.iterrows():
        name = str(r["target"])
        rec = {}
        for sp_csv, sp_key in [("NaCl", "nacl"), ("KCl", "kcl"),
                                ("H2O", "h2o"), ("RRL", "rrl"),
                                ("SiO", None), ("SO", None), ("COM", "com")]:
            if sp_key is None:
                continue
            k = f"{sp_csv}_kind"
            if k in df.columns:
                v = str(r.get(k, "") or "").strip().lower()
                if v in ("det", "tent", "ul", "na"):
                    rec[sp_key] = v
        out[name] = rec
    return out


def detections_for_target(name):
    """Return dict {species: dict(observed=bool, detected=bool)}.

    None if no pipeline analysis at all.
    """
    adirs = sorted((ANALYSIS / name).glob("2*")) if (ANALYSIS / name).is_dir() else []
    if not adirs:
        return None
    obs = {sp: False for sp in SPECIES}
    det = {sp: False for sp in SPECIES}
    seen_analysis = False
    for ad in adirs:
        meas = ad / "line_measurements.csv"
        if not meas.exists():
            continue
        try:
            m = pd.read_csv(meas)
        except pd.errors.EmptyDataError:
            continue
        if m.empty or not {"source", "snr", "line"} <= set(m.columns):
            continue
        bid = brightest_id(ad)
        if bid is None:
            continue
        seen_analysis = True
        # Coverage: any row at any source with a known species
        for _, row in m.iterrows():
            sp = species_of(row.get("line"), row.get("group"))
            if sp:
                obs[sp] = True
        # Detection: row at brightest source with snr>=5
        sub = m[(m["source"] == bid) & (m["snr"] >= 5.0)]
        for _, row in sub.iterrows():
            sp = species_of(row.get("line"), row.get("group"))
            if sp:
                det[sp] = True
    if not seen_analysis:
        return None
    return {sp: dict(observed=obs[sp], detected=det[sp]) for sp in SPECIES}


def merge_lit(state, lit_rec):
    """Apply lit overrides: det → counts as detected+observed; ul → observed,
    not detected; na → not observed (treat as no coverage); tent → observed
    but NOT counted as detected (conservative)."""
    if not lit_rec:
        return state
    for sp, val in lit_rec.items():
        if sp not in state:
            continue
        if val == "det":
            state[sp]["observed"] = True
            state[sp]["detected"] = True
        elif val == "tent":
            state[sp]["observed"] = True
        elif val == "ul":
            state[sp]["observed"] = True
            # don't flip detected
        elif val == "na":
            state[sp]["observed"] = False
            state[sp]["detected"] = False
    return state


def parse_lt3kpc(md_path):
    names = []
    in_table = False
    with md_path.open() as fh:
        for line in fh:
            if line.startswith("| # |"):
                in_table = True
                continue
            if not in_table:
                continue
            m = re.match(r"^\|\s*\d+\s*\|\s*([^|]+)\|", line)
            if m:
                names.append(m.group(1).strip())
    return names


def summarize(sample_df, label, lit):
    n = len(sample_df)
    n_alma = int((sample_df["n_alma_obs"].fillna(0) > 0).sum())
    analyzed = 0
    n_obs = {sp: 0 for sp in SPECIES}
    n_det = {sp: 0 for sp in SPECIES}
    for _, r in sample_df.iterrows():
        d = detections_for_target(r["name"])
        if d is None:
            # Allow lit-only entries even without pipeline (e.g., Orion-SrcI)
            lit_rec = lit.get(r["name"])
            if lit_rec:
                d = {sp: dict(observed=False, detected=False) for sp in SPECIES}
            else:
                continue
        d = merge_lit(d, lit.get(r["name"]))
        analyzed += 1
        for sp in SPECIES:
            if d[sp]["observed"]:
                n_obs[sp] += 1
            if d[sp]["detected"]:
                n_det[sp] += 1
    return dict(label=label, n=n, n_alma=n_alma, analyzed=analyzed,
                n_obs=n_obs, n_det=n_det)


def cell(n_det, n_obs):
    if n_obs == 0:
        return r"\nodata"
    pct = 100.0 * n_det / n_obs
    return f"{pct:.0f}\\% ({n_det}/{n_obs})"


def main():
    src = pd.read_csv(L4_CSV)
    lit = load_lit()
    s2 = src[src["dist_kpc"] <= 2.0].copy()
    print(f"2 kpc sample (from L4_d2 csv): {len(s2)}")

    names_3 = parse_lt3kpc(MD3)
    print(f"3 kpc sample from MD: {len(names_3)} names")
    have_l4 = {r["name"]: r for _, r in src.iterrows()}
    s3_rows = []
    for nm in names_3:
        if nm in have_l4:
            s3_rows.append(have_l4[nm])
        else:
            s3_rows.append(pd.Series(dict(name=nm, dist_kpc=np.nan,
                                           lbol_lsun=np.nan, n_alma_obs=0)))
    s3 = pd.DataFrame(s3_rows)
    print(f"3 kpc sample combined: {len(s3)}")

    rows = [summarize(s2, "$d \\leq 2$\\,kpc", lit),
            summarize(s3, "$d \\leq 3$\\,kpc", lit)]
    for r in rows:
        print(r)

    out = []
    out.append(r"\begin{deluxetable}{lcccccccccc}")
    out.append(r"\tabletypesize{\footnotesize}")
    out.append(r"\tablecaption{Demographic summary of the L$_\mathrm{bol} \geq 10^4\,L_\odot$ "
               r"HMYSO sample within $d \leq 2$\,kpc and $d \leq 3$\,kpc. "
               r"`Analyzed' counts targets with a completed "
               r"\texttt{line\_pipeline} run or a literature-confirmed "
               r"detection. Per-species cells give "
               r"$N_\mathrm{det}/N_\mathrm{obs}$ as a percentage, where "
               r"$N_\mathrm{obs}$ is the number of analyzed targets whose "
               r"spectral coverage included that species and "
               r"$N_\mathrm{det}$ is the subset with a $\geq 5\sigma$ "
               r"detection at the brightest mm continuum source or a "
               r"literature-confirmed detection. The denominator differs "
               r"between species because not every ALMA program covered "
               r"every line.\label{tab:samplesummary}}")
    out.append(r"\tablehead{")
    out.append(r"\colhead{Sample} & \colhead{$N$} & "
               r"\colhead{$N_\mathrm{ALMA}$} & \colhead{Analyzed} & "
               r"\colhead{RRL} & \colhead{NaCl} & \colhead{KCl} & "
               r"\colhead{H$_2$O\,232} & \colhead{COM} & \colhead{SiS} }")
    out.append(r"\startdata")
    for r in rows:
        # n_alma + analyzed cells: pct of N
        n_alma_cell = (f"{r['n_alma']} "
                       f"({100*r['n_alma']/max(r['n'],1):.0f}\\%)")
        an_cell = (f"{r['analyzed']} "
                   f"({100*r['analyzed']/max(r['n'],1):.0f}\\%)")
        out.append(
            f"{r['label']} & {r['n']} & {n_alma_cell} & {an_cell} & "
            f"{cell(r['n_det']['rrl'], r['n_obs']['rrl'])} & "
            f"{cell(r['n_det']['nacl'], r['n_obs']['nacl'])} & "
            f"{cell(r['n_det']['kcl'], r['n_obs']['kcl'])} & "
            f"{cell(r['n_det']['h2o'], r['n_obs']['h2o'])} & "
            f"{cell(r['n_det']['com'], r['n_obs']['com'])} & "
            f"{cell(r['n_det']['sis'], r['n_obs']['sis'])} \\\\"
        )
    out.append(r"\enddata")
    out.append(r"\tablecomments{$N_\mathrm{ALMA}$ counts targets with at "
               r"least one ALMA observation in the archive whose angular "
               r"resolution corresponds to $<500$\,AU at the source "
               r"distance. Percentages in the $N_\mathrm{ALMA}$ and "
               r"Analyzed columns are with respect to $N$. Per-species "
               r"percentages are with respect to "
               r"$N_\mathrm{obs}$ for that species. Literature detections "
               r"are taken from data/literature\_detections.csv "
               r"(Ginsburg+2019 Orion-SrcI; Plambeck+2013 Orion-BN; "
               r"and future entries).}")
    out.append(r"\end{deluxetable}")
    OUT_TEX.write_text("\n".join(out) + "\n")
    print(f"\nwrote {OUT_TEX}")


if __name__ == "__main__":
    main()
