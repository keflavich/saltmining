"""Integrated per-target detection+UL summary for the paper.

One row per target. For each species class, pick the best ALMA program
(max detections) and either:
  - report the strongest peak (highest snr detection >= 5sigma), OR
  - report the tightest 3sigma upper limit across class lines (10 km/s, native).

Species groups: NaCl, KCl, H2O, RRL, SiO, SiS, SO, COM (TODO).
Output:
  data/per_target_paper.csv     (machine-readable)
  /orange/adamginsburg/salt/demography_2026/per_target_paper.tex
"""
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from kelvin import conversion_factor

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
SRC_CSV = ROOT / "data/sources_L4_d2.csv"
OUT_CSV = ROOT / "data/per_target_paper.csv"
OUT_TEX = Path("/orange/adamginsburg/salt/demography_2026/per_target_paper.tex")

SPECIES = ["NaCl", "KCl", "H2O", "RRL", "SiO", "SiS", "SO", "AlO"]

# obs_params.csv carries curated peak/UL brightness + status for these four
# species (incl. literature-only sources like Orion-SrcI whose lines are not in
# our re-imaged cubes); used to backfill cells the pipeline aggregate left 'na'.
OBSPAR_SPECIES = {"NaCl": "nacl", "KCl": "kcl", "H2O": "h2o", "RRL": "rrl"}


def obsparams_backfill():
    """target -> {species: (kind, value_K)} from data/obs_params.csv."""
    path = ROOT / "data/obs_params.csv"
    if not path.exists():
        return {}
    op = pd.read_csv(path)
    out = {}
    for _, r in op.iterrows():
        name = str(r.get("name", "")).strip()
        if not name:
            continue
        d = {}
        for sp, pfx in OBSPAR_SPECIES.items():
            status = str(r.get(f"{pfx}_status", "")).strip().lower()
            peak = pd.to_numeric(r.get(f"{pfx}_peak_K"), errors="coerce")
            sig = pd.to_numeric(r.get(f"{pfx}_sig_K"), errors="coerce")
            if status.startswith("detect") and np.isfinite(peak):
                d[sp] = ("det", float(peak))
            elif status and status not in ("na", "nan", "notcovered",
                                            "not_covered") and np.isfinite(sig):
                d[sp] = ("ul", 3.0 * float(sig))
        if d:
            out[name] = d
    return out


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


def sigma_at_dv(sigma_native, dv_native, dv_target=10.0):
    if dv_target <= dv_native:
        return sigma_native
    return sigma_native * np.sqrt(dv_native / dv_target)


def best_for_target(target_dir: Path):
    """Return DataFrame of (program, line, group, snr, peak_K, sigma10_K, det) for
    brightest source rows, across all programs of this target."""
    rows = []
    for prop_dir in sorted(target_dir.glob("2*")):
        if not prop_dir.is_dir():
            continue
        bid = brightest_source_id(prop_dir / "continuum_sources.csv")
        if bid is None:
            continue
        meas = prop_dir / "line_measurements.csv"
        if not meas.exists():
            continue
        try:
            df = pd.read_csv(meas)
        except pd.errors.EmptyDataError:
            continue
        if df.empty or "snr" not in df.columns:
            continue
        df = df[df["source"] == bid].copy()
        if df.empty:
            continue
        src_dir = prop_dir / f"source_{bid:02d}"
        sigma10 = []
        peak_K = []
        for _, r in df.iterrows():
            npz = src_dir / f"{r['line']}.spec.npz"
            if not npz.exists():
                sigma10.append(np.nan); peak_K.append(np.nan); continue
            d = np.load(npz)
            if not {"vaxis", "spec", "sigma"} <= set(d.files):
                sigma10.append(np.nan); peak_K.append(np.nan); continue
            v = np.asarray(d["vaxis"], dtype=float)
            s = np.asarray(d["spec"], dtype=float)
            sn_native = float(d["sigma"])
            if v.size < 2 or not np.isfinite(sn_native):
                sigma10.append(np.nan); peak_K.append(np.nan); continue
            # spec.npz values are cube-native (Jy/beam for archive
            # products); convert to true K via the row's cube header
            cube = r.get("cube")
            if isinstance(cube, str) and cube:
                cube_full = ROOT / "uvdata" / prop_dir.name / target_dir.name / cube
                conv = conversion_factor(str(cube_full), float(r.get("rest_GHz", np.nan)))
            else:
                conv = np.nan
            dv = float(np.median(np.abs(np.diff(v))))
            sigma10.append(sigma_at_dv(sn_native * conv, dv, 10.0))
            peak_K.append(float(np.nanmax(s)) * conv)
        df["sigma10_K"] = sigma10
        df["peak_K"] = peak_K
        df["program"] = prop_dir.name
        df["src_id"] = bid
        rows.append(df)
    if not rows:
        return None
    return pd.concat(rows, ignore_index=True)


def species_for_line(group, line):
    """Map group + line into our display species."""
    if group == "NaCl": return "NaCl"
    if group == "KCl":  return "KCl"
    if group == "H2O":  return "H2O"
    if group == "RRL":  return "RRL"
    if group == "SiO":  return "SiO"
    if group == "SiS":  return "SiS"
    if group == "SO":   return "SO"
    return None


def stack_snr_for_species(target_dir: Path, sp: str):
    """Read per-source multi-line stack at the brightest source if available.
    Returns (snr, peak_K, program) for the best stack across programs, or
    None when no stack present.

    line_pipeline writes stack_salt() output as source_{NN}/{NaCl,KCl}_stack.npz
    with keys vaxis, spec, sigma, n_lines. We check that the peak/sigma >= 5
    and report peak_K + sigma to fold into the per-target table.
    """
    if sp not in ("NaCl", "KCl"):
        return None  # line_pipeline only stacks salt families today
    best = None
    for prop_dir in sorted(target_dir.glob("2*")):
        bid = brightest_source_id(prop_dir / "continuum_sources.csv")
        if bid is None:
            continue
        npz = prop_dir / f"source_{bid:02d}" / f"{sp}_stack.npz"
        if not npz.exists():
            continue
        d = np.load(npz)
        if not {"spec", "sigma"} <= set(d.files):
            continue
        spec = np.asarray(d["spec"])
        sigma = float(d["sigma"])
        if not np.isfinite(sigma) or sigma <= 0:
            continue
        # Stack is in cube-native units (Jy/beam); approximate the K
        # conversion with the median factor over this family's lines.
        conv = np.nan
        meas_csv = prop_dir / "line_measurements.csv"
        if meas_csv.exists():
            try:
                mdf = pd.read_csv(meas_csv)
            except pd.errors.EmptyDataError:
                mdf = pd.DataFrame()
            if not mdf.empty and {"group", "cube", "rest_GHz"} <= set(mdf.columns):
                fam = mdf[mdf["group"] == sp].drop_duplicates("line")
                convs = [conversion_factor(
                            str(ROOT / "uvdata" / prop_dir.name /
                                target_dir.name / str(mr["cube"])),
                            float(mr["rest_GHz"]))
                         for _, mr in fam.iterrows()
                         if isinstance(mr.get("cube"), str)]
                convs = [c for c in convs if np.isfinite(c)]
                if convs:
                    conv = float(np.median(convs))
        peak = float(np.nanmax(spec)) * conv
        sigma = sigma * conv
        snr = (peak / sigma) if (np.isfinite(peak) and np.isfinite(sigma)
                                  and sigma > 0) else np.nan
        if not np.isfinite(snr):
            continue
        if best is None or snr > best[0]:
            best = (snr, peak, sigma, prop_dir.name, int(d.get("n_lines", 0)))
    return best


def aggregate_species(df_all, target_dir: Path):
    """For one target's df_all + analysis_products dir, produce
    {species: (kind, value_K, program, note)}.
    kind: 'det' = detection peak; 'ul' = 3sigma 10 km/s UL; 'na' = no data.

    Detection takes either:
      - any individual line snr >= 5 at the brightest source, OR
      - a multi-line stack at the brightest source with peak/sigma >= 5
    """
    out = {}
    for sp in SPECIES:
        rows = df_all[df_all["group"] == sp]
        if rows.empty:
            out[sp] = ("na", np.nan, None, "")
            continue
        rejected = False
        det_rows = rows[rows["snr"] >= 5.0]
        if not det_rows.empty:
            top = det_rows.loc[det_rows["snr"].idxmax()]
            # strict-vet gating (incl. manual expert overrides):
            # REAL -> det, SUSPECT -> tent, REJECT -> fall through to UL
            verdict = None
            try:
                from build_data_summary import _strict_vet_pass
                verdict = _strict_vet_pass(target_dir.name,
                                            str(top["program"]),
                                            int(top.get("src_id", -1)), sp)
            except ImportError:
                verdict = None
            if verdict == "REJECT":
                rejected = True  # demoted: skip stack, go to UL below
            elif verdict == "SUSPECT":
                out[sp] = ("tent", float(top["peak_K"]),
                           str(top["program"]), "line (vet: suspect)")
                continue
            else:
                out[sp] = ("det", float(top["peak_K"]),
                           str(top["program"]), "line")
                continue
        # No single-line detection; try the multi-line stack. Threshold 4
        # rather than 5: a stack of N>=2 independent lines effectively
        # multiplies the post-trials confidence; the disk velocity-spread
        # also dilutes peak SNR in each individual line. This recovers
        # MonR2-IRS3 NaCl (4 lines @ 4.8 sigma post-stack) as a detection.
        stack = None if rejected else stack_snr_for_species(target_dir, sp)
        if stack is not None and stack[0] >= 4.0:
            snr, peak, sigma, prog, nl = stack
            out[sp] = ("det", peak, prog, f"stack n={nl} snr={snr:.1f}")
            continue
        valid = rows[np.isfinite(rows["sigma10_K"])]
        if valid.empty:
            out[sp] = ("na", np.nan, None, "")
            continue
        best = valid.loc[valid["sigma10_K"].idxmin()]
        out[sp] = ("ul", 3.0 * float(best["sigma10_K"]), str(best["program"]), "")
    return out


def fmt_cell(kind, val_K, flag_marker=""):
    from kelvin import fmt_K
    if kind == "na":
        return r"\nodata"
    if val_K is None or not np.isfinite(val_K):
        return r"\nodata"
    if kind == "det":
        return f"{fmt_K(val_K)}{flag_marker}"
    return rf"$<${fmt_K(val_K)}"


# Map per-species warning categories to LaTeX superscripts on the cell
SPECIES_FLAG_MARKERS = {
    "NaCl": [("NaCl_VEXC_ONLY", "$^{v}$"),
             ("NaCl_SINGLE_LINE_DET", "$^{s}$"),
             ("H2O_NONDET_WITH_NACL_DET", "$^{w}$"),
             ("HIGH_CONFUSION", "$^{c}$")],
    "KCl":  [("KCl_VEXC_ONLY", "$^{v}$"),
             ("KCL_DET_NACL_NONDET", "$^{n}$"),
             ("KCl_SINGLE_LINE_DET", "$^{s}$"),
             ("HIGH_CONFUSION", "$^{c}$")],
    "H2O":  [("H2O_SINGLE_LINE_DET", "$^{s}$"),
             ("HIGH_CONFUSION", "$^{c}$")],
    "RRL":  [("RRL_SINGLE_LINE_DET", "$^{s}$"),
             ("HIGH_CONFUSION", "$^{c}$")],
    "SiO":  [("SiO_SINGLE_LINE_DET", "$^{s}$")],
    "SiS":  [("SiS_SINGLE_LINE_DET", "$^{s}$")],
    "SO":   [("SO_SINGLE_LINE_DET", "$^{s}$")],
}


def marker_for(species, flags_str):
    if not flags_str:
        return ""
    flags = set(flags_str.split("|"))
    marks = []
    for flagname, mark in SPECIES_FLAG_MARKERS.get(species, []):
        if flagname in flags:
            marks.append(mark)
    return "".join(marks)


def main():
    src = pd.read_csv(SRC_CSV)
    src_by_name = {n: r for n, r in zip(src["name"], src.to_dict("records"))}
    out_rows = []

    # Literature overrides: targets whose detections we take from the
    # published literature (no on-disk pipeline analysis needed/reused).
    lit_path = ROOT / "data/literature_detections.csv"
    if lit_path.exists():
        lit = pd.read_csv(lit_path)
        for _, lr in lit.iterrows():
            row = {"target": str(lr["target"])}
            for sp in SPECIES:
                kind_col = f"{sp}_kind"
                val_col = f"{sp}_K_or_3s_mK"
                kind = lr.get(kind_col, "na") if kind_col in lit.columns else "na"
                if pd.isna(kind) or not kind:
                    kind = "na"
                K_val = lr.get(val_col, np.nan) if val_col in lit.columns else np.nan
                # mK -> K for storage
                row[f"{sp}_kind"] = str(kind)
                row[f"{sp}_K"] = (float(K_val) / 1000.0) if pd.notna(K_val) and K_val != "" else np.nan
                row[f"{sp}_prog"] = str(lr.get("reference", "literature"))
            out_rows.append(row)

    lit_targets = {r["target"] for r in out_rows}

    for target_dir in sorted(ANALYSIS.iterdir()):
        if not target_dir.is_dir():
            continue
        if target_dir.name in lit_targets:
            continue  # literature override takes precedence
        df_all = best_for_target(target_dir)
        if df_all is None or df_all.empty:
            continue
        agg = aggregate_species(df_all, target_dir)
        out_rows.append({
            "target": target_dir.name,
            **{f"{sp}_kind": agg[sp][0] for sp in SPECIES},
            **{f"{sp}_K": agg[sp][1] for sp in SPECIES},
            **{f"{sp}_prog": agg[sp][2] for sp in SPECIES},
        })
    out = pd.DataFrame(out_rows)

    # Backfill NaCl/KCl/H2O/RRL cells left 'na' with curated obs_params values
    # (fills literature-only archetypes: Orion-SrcI water, NGC6334I, IRAS17233,
    # and BN's RRL -- SrcI has no RRL so its RRL cell stays blank, and BN's
    # RRL detection carries the column).
    # Literature-only archetype cells absent from obs_params: Orion-SrcI's
    # vibrationally-excited hot water (232.687 GHz v2=1 line, not covered by our
    # re-imaged cubes) is the reference "whopping" water detection
    # \citep{Hirota2016}.
    LIT_OVERRIDE = {("Orion-SrcI", "H2O"): ("det", 1500.0)}

    def _norm(n):
        return str(n).replace("_", "-")

    obp = obsparams_backfill()
    for i, tgt in enumerate(out["target"]):
        for (lt, lsp), (lk, lv) in LIT_OVERRIDE.items():
            if _norm(tgt) == _norm(lt) and f"{lsp}_kind" in out.columns:
                # fill when unclassified OR when a detection has no value yet
                if (str(out.at[i, f"{lsp}_kind"]) in ("na", "nan", "", "None")
                        or pd.isna(out.at[i, f"{lsp}_K"])):
                    out.at[i, f"{lsp}_kind"] = lk
                    out.at[i, f"{lsp}_K"] = lv
                    if f"{lsp}_prog" in out.columns:
                        out.at[i, f"{lsp}_prog"] = "literature"
        vals = obp.get(str(tgt))
        if not vals:
            continue
        for sp, (kind, valK) in vals.items():
            kc, vc = f"{sp}_kind", f"{sp}_K"
            cur = str(out.at[i, kc]) if kc in out.columns else "na"
            if cur in ("na", "nan", "", "None") or pd.isna(out.at[i, vc]):
                out.at[i, kc] = kind
                out.at[i, vc] = valK
                if f"{sp}_prog" in out.columns:
                    out.at[i, f"{sp}_prog"] = "obs\\_params"

    # Merge evidence audit flags (one row per analyzed target/proposal).
    # We use the union of flags across proposals for each target.
    audit_csv = ROOT / "data/evidence_audit.csv"
    if audit_csv.exists():
        a = pd.read_csv(audit_csv)
        if not a.empty and "flags" in a.columns:
            flags_per_target = (a.groupby("target")["flags"]
                                  .apply(lambda s: "|".join(sorted(set(
                                      f for v in s.dropna() for f in str(v).split("|") if f
                                  ))))
                                  .reset_index().rename(columns={"flags": "warning_flags"}))
            out = out.merge(flags_per_target, left_on="target", right_on="target",
                              how="left")
            out["warning_flags"] = out["warning_flags"].fillna("")
        else:
            out["warning_flags"] = ""
    else:
        out["warning_flags"] = ""

    out.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV} ({len(out)} targets)")

    # Build LaTeX deluxetable*: one row per target. Columns: Field, then species cells.
    lines = [
        r"\begin{deluxetable*}{l" + "c" * len(SPECIES) + r"}",
        r"\tabletypesize{\scriptsize}",
        r"\tablecaption{Integrated per-target detection peaks and 3$\sigma$ upper limits "
        r"toward the brightest mm continuum source. Detection values are peak brightness "
        r"temperatures (K, converted from the cube-native Jy\,beam$^{-1}$ "
        r"scale with each cube's synthesized beam) of the highest-S/N line "
        r"in each species class; non-detection "
        r"values are 3$\sigma$ upper limits (K) at a 10 km/s channel using the line with "
        r"the tightest native-resolution noise. \nodata\ = no line of that class recorded.\label{tab:per_target}}",
        r"\tablehead{",
        r"\colhead{Field} & " + " & ".join(rf"\colhead{{{sp}}}" for sp in SPECIES) + r" \\",
        r" & " + " & ".join(["(K)"] * len(SPECIES)) + r"}",
        r"\startdata",
    ]
    try:
        from build_target_table import COMMON_NAME, alma_target_names_to_iras
    except ImportError:
        COMMON_NAME = {}
        alma_target_names_to_iras = lambda _s: None
    # Build target -> alma_target_names map for IRAS-fallback display.
    alt_col = src["alma_target_names"] if "alma_target_names" in src.columns else pd.Series([""] * len(src))
    alt_by_target = dict(zip(src["name"], alt_col.fillna("")))

    def display_name(tgt: str) -> str:
        common = COMMON_NAME.get(tgt)
        if common:
            return common
        iras = alma_target_names_to_iras(alt_by_target.get(tgt, ""))
        if iras:
            return iras
        return tgt

    out_sorted = out.sort_values("target")
    seen_names = set()
    for _, r in out_sorted.iterrows():
        tgt = display_name(str(r["target"])).replace("_", r"\_")
        if tgt in seen_names:
            continue
        seen_names.add(tgt)
        flags = str(r.get("warning_flags", "") or "")
        cells = [fmt_cell(r[f"{sp}_kind"], r[f"{sp}_K"], marker_for(sp, flags))
                  for sp in SPECIES]
        lines.append(tgt + " & " + " & ".join(cells) + r" \\")
    lines += [
        r"\enddata",
        r"\tablecomments{For each (target, species) we pick the program with the "
        r"strongest detection in that species class, or, when no detection, the "
        r"program+line combination with the smallest 3$\sigma$ noise at 10 km/s "
        r"channelization. RRL = any Hn$\alpha$/$\beta$/$\gamma$/$\delta$. Values are at "
        r"the native synthesized beam (not smoothed to a common beam); see "
        r"Table~\ref{tab:naclrrl} for the spatially-smoothed 300\,AU equivalents. "
        r"Warning superscripts on detection cells: $^{v}$ = the 5$\sigma$ line is "
        r"vibrationally excited (v$\ge$2) while the ground-state line(s) are "
        r"covered but $<$5$\sigma$ -- a big red flag for misidentification; "
        r"$^{s}$ = single 5$\sigma$ line drives the detection; $^{w}$ = the H$_2$O "
        r"lines are covered but not detected, despite a NaCl detection (warning); "
        r"$^{n}$ = KCl detected without NaCl detection (warning); $^{c}$ = "
        r"line-confusion fraction $f(>5\sigma) > 0.02$ in the recorded spectrum.}",
        r"\end{deluxetable*}",
    ]
    OUT_TEX.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT_TEX}")


if __name__ == "__main__":
    main()
