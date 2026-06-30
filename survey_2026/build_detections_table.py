"""Auto-generate the detections.tex summary table from per_target_paper.csv
+ literature_detections.csv.

Replaces the historically-hand-curated detections.tex (tab:similarities).

For each target, one row with yes/no/?/yes*/\nodata cells per species.
Species columns: disk, H$_2$O, NaCl, KCl, SiO, RRL, COMs, SiS, SO, PN.
- disk: from disk_classification.csv (manual override) — TODO/?
- H2O, NaCl, KCl, SiO, RRL, SiS, SO: from per_target_paper.csv kind column
- COMs, PN: pipeline cannot reliably classify yet — manual_overrides.csv

Cell translation:
   det  -> yes
   tent -> yes*
   ul   -> no
   na   -> \nodata
   ""   -> ?
   "?"  -> ?

Source label = literature mm/SMA name when known, else COMMON_NAME, else
target name with underscores escaped.
"""
import warnings
from pathlib import Path

import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
PER_TARGET = ROOT / "data/per_target_paper.csv"
LIT_CSV = ROOT / "data/literature_detections.csv"
DISK_CSV = ROOT / "data/disk_classification.csv"
OUT_TEX = Path("/orange/adamginsburg/salt/demography_2026/detections.tex")

SPECIES_COLS = ["H2O", "NaCl", "KCl", "SiO", "RRL", "SiS", "SO"]

# Override: literature COM and PN detections that the pipeline cannot
# classify. Keyed by target (CSV name).
MANUAL = {
    "Orion-SrcI":         dict(disk="yes",  COM="yes",  PN="?"),
    "Orion-BN":           dict(disk="no",   COM=r"\nodata", PN=r"\nodata"),
    "NGC6334I":           dict(disk="yes-c", COM="yes",  PN="yes*"),
    "NGC6334IN":          dict(disk="no",   COM="yes",  PN="yes*"),
    "MonR2-IRS3":         dict(disk="yes",  COM="?",    PN="?"),
    "MonR2-IRS2":         dict(disk="cont", COM="no",   PN="no"),
    "G336.4917-01.4741A": dict(disk="no",   COM="yes",  PN="yes*"),
    "G336.4917-01.4741B": dict(disk="yes*", COM="yes",  PN="yes*"),
    "G011.9197-00.6131":  dict(disk="yes",  COM="yes",  PN="yes*"),
    "G019.6097-00.2342":  dict(disk="yes-c", COM="yes", PN="no*"),
    "G023.0099-00.4108":  dict(disk="no",   COM="yes",  PN="no*"),
    "G017.6396+00.1580":  dict(disk="yes",  COM="no",   PN="no"),
    "G268.4222-00.8490":  dict(disk="?",    COM=r"\textsc{wip}", PN=r"\textsc{wip}"),
    "G326.6618+00.5207":  dict(disk="yes",  COM="yes*", PN="?"),
    "G345.5043+00.3480":  dict(disk="yes",  COM="yes",  PN="?"),
    "G345.4938+01.4677":  dict(disk="no",   COM="ext",  PN="no"),
    "G351.7745-00.5377":  dict(disk="yes-c", COM="no",  PN="yes"),
    "G013.6562-00.5997":  dict(disk="unres", COM="yes", PN="yes"),
    "G192.6005-00.0479":  dict(disk="no*",  COM="yes",  PN="no"),
    "G232.6207+00.9959":  dict(disk="yes",  COM="?",    PN="?"),
    "G015.0357-00.6795":  dict(disk="?",    COM="?",    PN="?"),
    "G189.0307+00.7821":  dict(disk="?",    COM="?",    PN="?"),
    "G345.0061+01.7944C": dict(disk="?",    COM="?",    PN="?"),
    "G345.0052+01.8209":  dict(disk="?",    COM="?",    PN="?"),
    "G010.8411-02.5919":  dict(disk="yes",  COM="yes",  PN="no"),
    "I16547-4247":        dict(disk="yes-c", COM="yes", PN="yes"),
}


def kind_to_cell(kind):
    if pd.isna(kind):
        return "?"
    k = str(kind).strip().lower()
    if k == "det":
        return "yes"
    if k == "tent":
        return "yes*"
    if k == "ul":
        return "no"
    if k == "na":
        return r"\nodata"
    return "?"


def _alt_names_lookup():
    """Build target -> alma_target_names map for IRAS-fallback display."""
    csv = ROOT / "data/sources_L4_d2.csv"
    if not csv.exists():
        return {}
    try:
        df = pd.read_csv(csv)
    except pd.errors.EmptyDataError:
        return {}
    alt = df["alma_target_names"] if "alma_target_names" in df.columns else pd.Series([""] * len(df))
    return dict(zip(df["name"], alt.fillna("")))


_ALT_BY_TARGET = None


def latex_source_name(name):
    """Display name: COMMON_NAME → IRAS XXXXX+YYYY → bare target."""
    global _ALT_BY_TARGET
    if _ALT_BY_TARGET is None:
        _ALT_BY_TARGET = _alt_names_lookup()
    try:
        from build_target_table import COMMON_NAME, alma_target_names_to_iras
    except ImportError:
        COMMON_NAME = {}
        alma_target_names_to_iras = lambda _s: None
    base = COMMON_NAME.get(name)
    if not base:
        base = alma_target_names_to_iras(_ALT_BY_TARGET.get(name, "")) or name
    return base.replace("_", r"\_")


def canonical_target(name):
    """Collapse alternate spellings: Orion_SrcI → Orion-SrcI etc."""
    aliases = {
        "Orion_SrcI": "Orion-SrcI",
        "OrionSrcI":  "Orion-SrcI",
        "OrionBN":    "Orion-BN",
    }
    return aliases.get(name, name)


def main():
    df = pd.read_csv(PER_TARGET)
    df["target"] = df["target"].apply(canonical_target)
    # Where literature and pipeline rows both exist, prefer the literature
    # version (kind in {det, tent, na}) over the pipeline (kind = ul).
    def is_lit(row):
        ks = [str(row.get(f"{sp}_kind", "")).lower() for sp in
              ["NaCl", "KCl", "H2O", "RRL", "SiO", "SiS", "SO"]]
        return any(k == "det" or k == "na" for k in ks) and "ul" not in ks
    df["__lit"] = df.apply(is_lit, axis=1)
    df = (df.sort_values(["target", "__lit"], ascending=[True, False])
          .drop_duplicates("target", keep="first")
          .drop(columns="__lit")
          .reset_index(drop=True))
    out = []
    out.append(r"\begin{deluxetable}{lcccccccccc}")
    out.append(r"\tabletypesize{\scriptsize}")
    out.append(r"\tablecaption{Summary of detections, tentative detections, "
               r"and non-detections of target species in the source sample. "
               r"Auto-generated from per\_target\_paper.csv (pipeline) + "
               r"literature\_detections.csv. yes = $\geq 5\sigma$ detection "
               r"or literature-confirmed detection at the brightest mm "
               r"continuum source; yes* = tentative; no = $\geq 5\sigma$ "
               r"upper limit; ? = ambiguous / not yet classified; "
               r"\nodata\ = species not covered in band.\label{tab:similarities}}")
    out.append(r"\tablehead{")
    out.append(r"\colhead{Source} & \colhead{disk} & \colhead{H$_2$O} & "
               r"\colhead{NaCl} & \colhead{KCl} & \colhead{SiO} & "
               r"\colhead{RRL} & \colhead{COMs} & \colhead{SiS} & "
               r"\colhead{SO} & \colhead{PN}}")
    out.append(r"\startdata")
    seen_names = set()
    for _, r in df.iterrows():
        tgt = str(r["target"])
        disp = latex_source_name(tgt)
        if disp in seen_names:
            continue
        seen_names.add(disp)
        manual = MANUAL.get(tgt, {})
        cells = [latex_source_name(tgt), manual.get("disk", "?")]
        for sp in SPECIES_COLS:
            cells.append(kind_to_cell(r.get(f"{sp}_kind")))
        cells.append(manual.get("COM", "?"))
        # Insert SiS+SO already done; reorder cells to match header
        # Header order: Source disk H2O NaCl KCl SiO RRL COMs SiS SO PN
        # Built order : Source disk H2O NaCl KCl SiO RRL SiS  SO  COMs
        # → swap COMs into position 8, then add SiS, SO, PN at 9, 10, 11
        # Simpler to rebuild in the right order:
        row_cells = {
            "disk": manual.get("disk", "?"),
            "H2O":  kind_to_cell(r.get("H2O_kind")),
            "NaCl": kind_to_cell(r.get("NaCl_kind")),
            "KCl":  kind_to_cell(r.get("KCl_kind")),
            "SiO":  kind_to_cell(r.get("SiO_kind")),
            "RRL":  kind_to_cell(r.get("RRL_kind")),
            "COMs": manual.get("COM", "?"),
            "SiS":  kind_to_cell(r.get("SiS_kind")),
            "SO":   kind_to_cell(r.get("SO_kind")),
            "PN":   manual.get("PN", "?"),
        }
        order = ["disk", "H2O", "NaCl", "KCl", "SiO", "RRL", "COMs",
                  "SiS", "SO", "PN"]
        line = " & ".join([latex_source_name(tgt)]
                          + [row_cells[k] for k in order]) + r" \\"
        out.append(line)
    out.append(r"\enddata")
    out.append(r"\tablecomments{`disk' = morphological or kinematic evidence "
               r"for a disk/continuum-source / unresolved (`yes-c'), `unres' "
               r"= unresolved, `cont' = continuum-only, `ext' = extended. "
               r"`yes' for any species means $\geq 5\sigma$ at the brightest "
               r"mm continuum source in at least one analyzed ALMA program, "
               r"OR a literature detection. Manual entries for `disk', "
               r"COMs, and PN columns are maintained in "
               r"build\_detections\_table.py until those modules can "
               r"auto-classify reliably.}")
    out.append(r"\end{deluxetable}")
    OUT_TEX.write_text("\n".join(out) + "\n")
    print(f"wrote {OUT_TEX} ({len(df)} rows)")


if __name__ == "__main__":
    main()
