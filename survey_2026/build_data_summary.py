"""Build the one-stop observation+detection summary table for demography_2026.

Columns (in order):
  Field             - IRAS-shorthand if any, else target name
  Source            - compact-source ID within the field
  theta_maj [AU]    - synthesized-beam major axis in AU at the source distance
  sigma_native      - per-channel noise at peak pixel of the brightest source (K)
  sigma_200au       - peak-T noise after spatial smoothing to 200 AU/dist arcsec (K)
  sigma_200au_10kms - sigma_200au at 10 km/s channel (K)
  f(>5sigma)        - line-confusion diagnostic (TODO)
  v_width           - cube velocity width per spectral window (km/s)
  d [kpc]           - adopted distance
  RRL               - any Hn alpha/beta/gamma >= 5 sigma in this field
  NaCl/H2O          - any NaCl, KCl, or H2O 232 line >= 5 sigma
  COM               - any complex-organic line (TODO)

Sources with no analysis_products row yet emit WIP/TODO placeholders.
"""
import json
import re
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
UVDIR = ROOT / "uvdata"
SRC_CSV = ROOT / "data/sources_L4_d2.csv"
SUMMARY_CSV = ROOT / "data/l4_d2_line_summary.csv"
INV_CSV = ROOT / "data/l4_d2_detection_inventory.csv"
OUT_TEX = Path("/orange/adamginsburg/salt/demography_2026/data_summary.tex")

WIP = r"\textsc{wip}"
TODO = r"\textsc{todo}"
NODATA = r"\nodata"

# Hand-curated IRAS shorthand for non-RMS or non-obvious cases.
IRAS_OVERRIDE = {
    "OrionBN-KL":         "I05327",
    "OrionB-Flame":       "I05393",
    "NGC6334I":           "I17175",
    "NGC6334IN":          "I17175N",
    "MonR2-IRS3":         "I06052",
    "MonR2-IRS2":         "I06052",
    "Lagoon-Her36":       "I18004",
    "NGC6514":            "I17590",
    "S140-IRS1":          "I22176",
    "GGD12-15":           "I06084",
    "IRAS17233-3606":     "I17233",
    "G353.2+0.9":         "I17220",
    "G189.0307+00.7821":  "I06056",
}

# Cube-name token -> IRAS shorthand, parsed from alma_target_names like
# "IRAS_07299-1651" -> "I07299".
# Matches "IRAS 09002-4732", "IRAS_09002-4732", or bare "I09002-4732".
IRAS_RE = re.compile(r"(?:IRAS[ _]*|\bI)([0-9]{5})(?:[+-][0-9]+)?", re.IGNORECASE)


def iras_shorthand(name: str, alt_names: str) -> str | None:
    """Return I-prefixed shorthand if available; None otherwise."""
    if name in IRAS_OVERRIDE:
        return IRAS_OVERRIDE[name]
    for tok in re.split(r"[;,]", str(alt_names or "")):
        m = IRAS_RE.search(tok)
        if m:
            return f"I{m.group(1)}"
    return None


def field_label(name: str, alt_names: str) -> str:
    """Common name with IRAS shorthand appended in parentheses when known."""
    iras = iras_shorthand(name, alt_names)
    if iras and iras != name:
        return f"{name}, {iras}"
    return name


# Per-target lit mm/SMA names keyed by (target, src_id).
# WARNING: src_id mapping varies by proposal; prefer CANONICAL_MM_COORDS
# (position-based) for multi-core sources.
LIT_MM_NAME = {
    ("I16547-4247", 1):  "I16547A",
    ("I16547-4247", 2):  "I16547B",
    ("Orion-SrcI", 1):   "Orion-SrcI",
    ("Orion-BN", 1):     "Orion-BN",
    ("OrionBN", 1):      "Orion-BN",
    ("MonR2-IRS3", 1):   "MonR2-IRS3",
    ("MonR2-IRS2", 1):   "MonR2-IRS2",
}

# Position-based canonical mm/SMA labels: target -> [(label, ra_deg, dec_deg, tol_arcsec)].
# When the brightest mm source in the analyzed proposal falls within tol of
# one of these positions, the source is labeled with that literature name.
CANONICAL_MM_COORDS = {
    "NGC6334I": [
        ("NGC6334I-MM1B", 260.22257, -35.78274, 0.5),
        ("NGC6334I-MM1D", 260.22270, -35.78284, 0.5),
        ("NGC6334I-MM1A", 260.22264, -35.78250, 0.5),
        ("NGC6334I-MM1C", 260.22260, -35.78262, 0.5),
        ("NGC6334I-MM2",  260.22126, -35.78427, 1.0),
        ("NGC6334I-MM3",  260.22243, -35.75109, 1.0),
    ],
    "NGC6334IN": [
        ("NGC6334I(N)-SMA6", 260.21796, -35.75476, 1.0),
        ("NGC6334I(N)-SMA1", 260.21804, -35.75517, 1.0),
        ("NGC6334I(N)-SMA2", 260.21788, -35.75539, 1.0),
        ("NGC6334I(N)-SMA4", 260.21788, -35.75488, 1.0),
    ],
}


def canonical_mm_label_for(target: str, ra_deg: float, dec_deg: float):
    """Match a (target, ra, dec) to a CANONICAL_MM_COORDS entry; return label
    or None."""
    coords = CANONICAL_MM_COORDS.get(target)
    if not coords:
        return None
    for label, lit_ra, lit_dec, tol_arcsec in coords:
        dra = (ra_deg - lit_ra) * np.cos(np.radians(dec_deg)) * 3600.0
        ddec = (dec_deg - lit_dec) * 3600.0
        if (dra * dra + ddec * ddec) ** 0.5 <= tol_arcsec:
            return label
    return None


def short_handle(name: str, alt_names: str) -> str:
    """A short token usable as the mm-naming prefix.
    IRAS shorthand if known, else the bare target name without the leading
    G... galactic-coord designation, else the name itself."""
    iras = iras_shorthand(name, alt_names)
    if iras and iras != name:
        return iras
    return name


_LIT_MM_SUFFIX_RE = re.compile(r"(?i)(?:-mm[0-9a-z]+|\(N\)-SMA[0-9]+|-IRS[0-9]+|-BN|-SrcI|/IRS[0-9]+|-A|-B|-C)$")


def source_label(target: str, src_id: int, alt_names: str,
                  brightness_rank: int) -> str:
    """Merged Source identifier for the data-summary table.
    Priority: literature mm name → common name (+ mm{rank} if name doesn't
    already encode a mm/component suffix) → IRASmm{rank} → targetmm{rank}.
    """
    key = (target, int(src_id))
    if key in LIT_MM_NAME:
        return LIT_MM_NAME[key]
    try:
        from build_target_table import COMMON_NAME
    except ImportError:
        COMMON_NAME = {}
    common = COMMON_NAME.get(target)
    if common:
        if _LIT_MM_SUFFIX_RE.search(common):
            return common if brightness_rank == 1 else f"{common}-mm{brightness_rank}"
        return f"{common}-mm{brightness_rank}"
    handle = short_handle(target, alt_names)
    return f"{handle}mm{brightness_rank}"


def brightness_rank_map(cont_csv: Path):
    """Return dict {src_id: rank} where rank 1 = brightest peak."""
    if not cont_csv.exists() or cont_csv.stat().st_size == 0:
        return {}
    try:
        df = pd.read_csv(cont_csv)
    except pd.errors.EmptyDataError:
        return {}
    if df.empty or "peak_Jybeam" not in df.columns:
        return {}
    ordered = df.sort_values("peak_Jybeam", ascending=False)["id"].tolist()
    return {int(sid): i + 1 for i, sid in enumerate(ordered)}


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


def beam_for_proposal(proposal_dir: Path):
    """Estimate representative cube beam in arcsec from the first cube
    referenced in line_measurements.csv."""
    csv = proposal_dir / "line_measurements.csv"
    if not csv.exists():
        return None
    try:
        df = pd.read_csv(csv)
    except pd.errors.EmptyDataError:
        return None
    if "cube" not in df.columns or df.empty:
        return None
    target = proposal_dir.parent.name
    proposal = proposal_dir.name
    for cube_name in df["cube"].dropna().unique():
        p = UVDIR / proposal / target / cube_name
        if not p.exists():
            continue
        try:
            h = fits.getheader(p, ext=0)
        except (OSError, fits.VerifyError):
            continue
        bmaj = h.get("BMAJ")
        if bmaj is not None:
            return float(bmaj) * 3600.0
        # Per-channel CASA cubes store beams in a BEAMS binary table
        # extension instead of the primary header.
        try:
            with fits.open(p, memmap=True) as hdul:
                for hdu in hdul[1:]:
                    if hdu.name == "BEAMS" and hdu.data is not None:
                        col = hdu.data["BMAJ"]
                        return float(np.nanmedian(col))  # arcsec
        except (OSError, fits.VerifyError, KeyError):
            continue
    return None


def native_noise_for(proposal_dir: Path, bid: int):
    """Median sigma_native (K) across all spec.npz under source_{bid}/."""
    sd = proposal_dir / f"source_{bid:02d}"
    if not sd.is_dir():
        return None, None
    sigs = []
    dvs = []
    for npz in sd.glob("*.spec.npz"):
        d = np.load(npz)
        if not {"vaxis", "sigma"} <= set(d.files):
            continue
        s = float(d["sigma"])
        v = np.asarray(d["vaxis"], dtype=float)
        if v.size < 2 or not np.isfinite(s):
            continue
        dv = float(np.median(np.abs(np.diff(v))))
        sigs.append(s); dvs.append(dv)
    if not sigs:
        return None, None
    return float(np.median(sigs)), float(np.median(dvs))


C_KMS = 2.99792458e5


def vwidth_for(proposal_dir: Path, bid: int):
    """Mean per-spw velocity span of the actual cubes used in this analysis
    (km/s). Reads NAXIS3 * CDELT3 / CRVAL3 * c from each cube header
    referenced by line_measurements.csv (the spec.npz axes only cover the
    pipeline's narrow cutout window and would underreport).
    """
    csv = proposal_dir / "line_measurements.csv"
    if not csv.exists():
        return None
    try:
        df = pd.read_csv(csv)
    except pd.errors.EmptyDataError:
        return None
    if "cube" not in df.columns or df.empty:
        return None
    target = proposal_dir.parent.name
    proposal = proposal_dir.name
    spans = []
    for cube_name in df["cube"].dropna().unique():
        cube_path = UVDIR / proposal / target / cube_name
        if not cube_path.exists():
            continue
        try:
            h = fits.getheader(cube_path, ext=0)
        except (OSError, fits.VerifyError):
            continue
        nz = h.get("NAXIS3")
        d_ = h.get("CDELT3")
        f0 = h.get("CRVAL3")
        if not (nz and d_ and f0):
            continue
        span_hz = abs(nz * d_)
        span_kms = span_hz / float(f0) * C_KMS
        spans.append(span_kms)
    return float(np.mean(spans)) if spans else None


def detections_at_brightest(proposal_dir: Path, bid: int):
    """Return dict with rrl, nacl, kcl, h2o, com booleans evaluated AT the
    brightest mm continuum source (snr >= 5 at source==bid). Returns None if
    measurements are missing.
    """
    csv = proposal_dir / "line_measurements.csv"
    if not csv.exists():
        return None
    try:
        df = pd.read_csv(csv)
    except pd.errors.EmptyDataError:
        return None
    if df.empty or not {"source", "line", "group", "snr"} <= set(df.columns):
        return None
    sub = df[(df["source"] == bid) & (df["snr"] >= 5.0)]
    rrl = bool(sub[(sub["group"] == "RRL")
                   | sub["line"].astype(str).str.match(r"^H\d+(alpha|beta|gamma|delta)$")
                   ].shape[0])
    nacl = bool(sub[sub["group"] == "NaCl"].shape[0])
    kcl = bool(sub[sub["group"] == "KCl"].shape[0])
    h2o = bool(sub[sub["group"] == "H2O"].shape[0])
    com_pat = re.compile(r"^(CH3OH|CH3CN|HC3N|CH3OCHO|C2H5CN|CH3OCH3|NH2CHO|HCOOH|HNCO)",
                         re.IGNORECASE)
    com = bool(sub[sub["line"].astype(str).apply(lambda L: bool(com_pat.match(L)))].shape[0])
    return dict(rrl=rrl, nacl=nacl, kcl=kcl, h2o=h2o, com=com)


SPECIES_LIT = ["NaCl", "KCl", "H2O", "SiO", "SiS", "SO", "COM", "RRL"]


def load_lit_refs():
    """Return {target: {ref: str, species_kind: {NaCl: det/ul/tent/na, ...}}}.
    Only targets with at least one 'det' or 'tent' species are returned.
    Adds underscore-/hyphen-swapped aliases so name spellings in sources CSV
    that differ from lit CSV still resolve."""
    lit_path = ROOT / "data/literature_detections.csv"
    if not lit_path.exists():
        return {}
    try:
        df = pd.read_csv(lit_path)
    except pd.errors.EmptyDataError:
        return {}
    out = {}
    for _, r in df.iterrows():
        target = str(r["target"])
        kinds = {s: str(r.get(f"{s}_kind", "") or "").strip().lower()
                  for s in SPECIES_LIT}
        if any(v in ("det", "tent") for v in kinds.values()):
            entry = dict(ref=str(r["reference"]), species_kind=kinds)
            out[target] = entry
            for alias in (target.replace("-", "_"), target.replace("_", "-")):
                if alias != target:
                    out.setdefault(alias, entry)
    return out


def apply_lit_override(species_key, mine_bool, lit_kinds):
    """Merge pipeline boolean with literature kind into a cell string.
    species_key: lit-side species name (NaCl/KCl/H2O/COM).
    mine_bool : True if our pipeline detects this species at bid.
    lit_kinds : dict from load_lit_refs entry, or {}.
    Returns one of: 'y', 'y?', 'n', '\\nodata'.
    """
    lit_kind = lit_kinds.get(species_key, "") if lit_kinds else ""
    if mine_bool:
        return "y"
    if lit_kind == "det":
        return r"y$^L$"
    if lit_kind == "tent":
        return r"y?$^L$"
    if lit_kind == "ul":
        return "n"
    if lit_kind == "na":
        return r"\nodata"
    return "n"


def _continuum_row(cont_csv: Path, src_id: int):
    if not cont_csv.exists() or cont_csv.stat().st_size == 0:
        return None
    try:
        df = pd.read_csv(cont_csv)
    except pd.errors.EmptyDataError:
        return None
    hit = df[df["id"] == src_id]
    if hit.empty:
        return None
    return hit.iloc[0]


def collect():
    src = pd.read_csv(SRC_CSV)
    lit_refs = load_lit_refs()

    # collect per-target best (target, proposal) by analysis_products
    rows = []
    seen_labels = set()
    for _, r in src.iterrows():
        name = r["name"]
        dist = r["dist_kpc"]
        merged_src_pre = source_label(name, 1, r.get("alma_target_names", ""), 1)
        # any analysis_products subdir for this target?
        adirs = sorted((ANALYSIS / name).glob("2*")) if (ANALYSIS / name).is_dir() else []
        # Pick the proposal with the most lines detected at its own brightest
        # mm continuum source; tie-break by smallest beam.
        best = None
        best_key = (-1, np.inf)  # (n_detections [neg later], theta_au)
        for ad in adirs:
            bid = brightest_source_id(ad / "continuum_sources.csv")
            if bid is None:
                continue
            beam_arcsec = beam_for_proposal(ad)
            theta_au = beam_arcsec * dist * 1000.0 if beam_arcsec else np.inf
            # Skip proposals whose resolution doesn't meet the <500 AU bar.
            if theta_au and np.isfinite(theta_au) and theta_au > 500.0:
                continue
            det = detections_at_brightest(ad, bid)
            n_det = 0
            if det is not None:
                n_det = sum(int(det[k]) for k in ("rrl", "nacl", "kcl", "h2o", "com"))
            # maximize n_det, then minimize theta_au
            key = (n_det, -theta_au)
            best_key_val = (best_key[0], -best_key[1])
            if key > best_key_val:
                best_key = (n_det, theta_au)
                best = (ad, bid, beam_arcsec, theta_au)
        if best is None:
            # No analysis_products in the expected location. Try alias
            # lookup (Orion-SrcI ↔ Orion_SrcI etc.) before giving up.
            for alias in (name.replace("-", "_"), name.replace("_", "-")):
                if alias == name:
                    continue
                adirs = sorted((ANALYSIS / alias).glob("2*")) if (ANALYSIS / alias).is_dir() else []
                if adirs:
                    name = alias  # use alias going forward
                    for ad in adirs:
                        bid = brightest_source_id(ad / "continuum_sources.csv")
                        if bid is None:
                            continue
                        beam_arcsec = beam_for_proposal(ad)
                        theta_au = beam_arcsec * dist * 1000.0 if beam_arcsec else np.inf
                        if theta_au and np.isfinite(theta_au) and theta_au > 500.0:
                            continue
                        det = detections_at_brightest(ad, bid)
                        n_det = sum(int(det[k]) for k in ("rrl","nacl","kcl","h2o","com")) if det else 0
                        if (n_det, -theta_au) > (best_key[0], -best_key[1]):
                            best_key = (n_det, theta_au)
                            best = (ad, bid, beam_arcsec, theta_au)
                    break
        if best is None:
            # Still nothing on disk. If lit reports a detection for this target
            # AND the source meets the <500 AU resolution requirement, emit a
            # lit-only WIP row. Else skip.
            lit_kinds = lit_refs.get(name, {}).get("species_kind", {})
            has_lit_det = any(v in ("det", "tent") for v in lit_kinds.values())
            br = r.get("best_res_arcsec")
            has_500au = pd.notna(br) and br * dist * 1000.0 < 500.0
            if not (has_lit_det or has_500au):
                continue
            if merged_src_pre in seen_labels:
                continue
            seen_labels.add(merged_src_pre)
            rows.append(dict(
                target=name,
                source=merged_src_pre,
                theta_au=TODO, sigma_native=TODO, sigma_200au=TODO,
                sigma_200au_10kms=TODO, f5sig=TODO,
                d=f"{dist:.2f}",
                rrl=WIP, nacl=WIP, kcl=WIP, h2o=WIP, com=WIP,
                lit_kinds=lit_kinds,
            ))
            continue
        ad, bid, beam_arcsec, theta_au = best
        # Position-based canonical mm label override (NGC6334I etc.)
        crow = _continuum_row(ad / "continuum_sources.csv", bid)
        canon_label = None
        if crow is not None and "ra_deg" in crow.index:
            canon_label = canonical_mm_label_for(
                name, float(crow["ra_deg"]), float(crow["dec_deg"]))
        rmap = brightness_rank_map(ad / "continuum_sources.csv")
        rank = rmap.get(int(bid), 1)
        merged_src = canon_label or source_label(
            name, bid, r.get("alma_target_names", ""), rank)
        # Dedupe: skip if label already emitted (e.g. two distinct G... targets
        # mapped to the same shorthand). The first occurrence wins.
        if merged_src in seen_labels:
            continue
        seen_labels.add(merged_src)
        sigma_nat, dv_nat = native_noise_for(ad, bid)
        if theta_au and np.isfinite(theta_au) and sigma_nat:
            if theta_au < 200.0:
                sigma_200 = sigma_nat * (theta_au / 200.0) ** 2
            else:
                sigma_200 = sigma_nat
            if dv_nat and dv_nat < 10.0:
                sigma_200_10 = sigma_200 * np.sqrt(dv_nat / 10.0)
            else:
                sigma_200_10 = sigma_200
        else:
            sigma_200 = None
            sigma_200_10 = None
        det = detections_at_brightest(ad, bid)
        if det is None:
            rrl = nacl = kcl = h2o = com = WIP
        else:
            rrl = "y" if det["rrl"] else "n"
            nacl = "y" if det["nacl"] else "n"
            kcl = "y" if det["kcl"] else "n"
            h2o = "y" if det["h2o"] else "n"
            com = "y" if det["com"] else "n"
        # aux: confusion fraction
        f5_str = TODO
        if (ROOT / "data/data_summary_aux.csv").exists():
            try:
                aux = pd.read_csv(ROOT / "data/data_summary_aux.csv")
                hit = aux[(aux["target"] == name) & (aux["proposal"] == ad.name)]
                if not hit.empty and pd.notna(hit.iloc[0]["f5sigma"]):
                    f5_str = f"{hit.iloc[0]['f5sigma']:.2f}"
                    if hit.iloc[0].get("com") is True and com == "n":
                        com = "y?"
            except pd.errors.EmptyDataError:
                pass
        lit_kinds = lit_refs.get(name, {}).get("species_kind", {})
        rows.append(dict(
            target=name,
            source=merged_src,
            theta_au=f"{theta_au:.0f}" if (theta_au and np.isfinite(theta_au)) else TODO,
            sigma_native=f"{sigma_nat*1000.0:.1f}" if sigma_nat else TODO,
            sigma_200au=f"{sigma_200*1000.0:.1f}" if sigma_200 else TODO,
            sigma_200au_10kms=f"{sigma_200_10*1000.0:.2f}" if sigma_200_10 else TODO,
            f5sig=f5_str,
            d=f"{dist:.2f}",
            rrl=rrl, nacl=nacl, kcl=kcl, h2o=h2o, com=com,
            lit_kinds=lit_kinds,
        ))
    return pd.DataFrame(rows)


def _latex_escape(s: str) -> str:
    return s.replace("&", r"\&").replace("_", r"\_").replace("%", r"\%")


def _classify_footnote(row, lit_refs):
    """Return (mark, ref_text or None).
    mark = 'a','b',... assigned via fn_map below if target has a lit detection;
    mark = '*' for new-this-work (any pipeline detection AND no lit-det);
    mark = None for non-detection / no information.
    """
    tgt = row["target"]
    lit_kinds = row.get("lit_kinds", {}) or {}
    has_lit_det = any(v in ("det", "tent") for v in lit_kinds.values())
    pipeline_det = any(row.get(k) == "y" for k in ("nacl", "kcl", "h2o", "com"))
    if has_lit_det:
        return ("lit", lit_refs.get(tgt, {}).get("ref", ""))
    if pipeline_det:
        return ("new", None)
    return (None, None)


def _merge_cell(species_key, mine, lit_kinds):
    """Merge per-species mine ('y'/'n'/WIP/TODO) with lit kind into final cell."""
    if mine in (WIP, TODO):
        return mine
    lit_kind = (lit_kinds or {}).get(species_key, "")
    if mine == "y":
        return "y"
    if lit_kind == "det":
        return r"y$^L$"
    if lit_kind == "tent":
        return r"y?$^L$"
    if lit_kind == "na":
        return r"\nodata"
    # default ('ul' or absent) -> show pipeline value
    return mine


def write_tex(df: pd.DataFrame):
    lit_refs = load_lit_refs()
    # Assign footnote letters to lit-references (in order of first appearance).
    fn_letter = {}
    star_used = False
    out = []
    out.append(r"\startlongtable")
    out.append(r"\begin{deluxetable}{lccccccccccc}")
    out.append(r"\tabletypesize{\scriptsize}")
    out.append(r"\tablecaption{Observation summary and detection inventory across the "
               r"full target sample. Beam in AU at the source distance, $\sigma$ in mK. "
               r"Detection columns: y = $\geq 5\sigma$ in our pipeline; "
               r"y$^L$ = literature detection (no pipeline coverage / below "
               r"$5\sigma$ here); n = non-detection; "
               r"\textsc{wip}\ = analysis in progress; \nodata\ = species not in band. "
               r"Source identifier is the literature mm/SMA name (matched by "
               r"position) when known, else \texttt{<common-name>-mm<rank>}, "
               r"ranking by mm peak brightness within the field. "
               r"Source-name footnotes link to first-detection references; "
               r"$^*$ marks first detection reported in this work.\label{tab:obs}}")
    out.append(r"\tablehead{")
    out.append(r"\colhead{Source} & "
               r"\colhead{$\theta_\mathrm{maj}$} & "
               r"\colhead{$\sigma_\mathrm{nat}$} & "
               r"\colhead{$\sigma_{200\mathrm{au}}$} & "
               r"\colhead{$\sigma_{200,10}$} & "
               r"\colhead{$f({>}5\sigma)$} & "
               r"\colhead{$d$} & "
               r"\colhead{RRL} & \colhead{NaCl} & \colhead{KCl} & "
               r"\colhead{H$_2$O} & \colhead{COM} \\")
    out.append(r" & (AU) & (mK) & (mK) & (mK) & & (kpc) & & & & 232 & }")
    out.append(r"\startdata")
    for _, r in df.iterrows():
        kind, ref = _classify_footnote(r, lit_refs)
        mark = ""
        if kind == "lit":
            tgt = r["target"]
            if tgt not in fn_letter:
                fn_letter[tgt] = chr(ord('a') + len(fn_letter))
            mark = f"\\tablenotemark{{{fn_letter[tgt]}}}"
        elif kind == "new":
            mark = r"\tablenotemark{*}"
            star_used = True
        src_disp = _latex_escape(str(r["source"])) + mark
        lit_k = r.get("lit_kinds", {}) or {}
        nacl = _merge_cell("NaCl", r["nacl"], lit_k)
        kcl  = _merge_cell("KCl",  r["kcl"],  lit_k)
        h2o  = _merge_cell("H2O",  r["h2o"],  lit_k)
        com  = _merge_cell("COM",  r["com"],  lit_k)
        out.append(
            f"{src_disp} & "
            f"{r['theta_au']} & {r['sigma_native']} & "
            f"{r['sigma_200au']} & {r['sigma_200au_10kms']} & "
            f"{r['f5sig']} & {r['d']} & "
            f"{r['rrl']} & {nacl} & {kcl} & {h2o} & {com} \\\\"
        )
    out.append(r"\enddata")
    out.append(r"\tablecomments{$\sigma_\mathrm{nat}$ is the per-channel noise at the "
               r"brightest mm continuum source's peak pixel. $\sigma_{200\mathrm{au}}$ "
               r"is the point-source-scaled noise at a $200$\,AU effective beam at the "
               r"source distance: $\sigma_T(200) = \sigma_T(\mathrm{nat}) \times "
               r"(\theta_\mathrm{nat}/\theta_{200})^2$ when $\theta_\mathrm{nat} < 200$\,AU; "
               r"otherwise the native value is reported. $\sigma_{200,10}$ is the same "
               r"after smoothing to a $10$\,\kms\ channel. $f({>}5\sigma)$ is the "
               r"fraction of channels above five times the field's robust noise; a high "
               r"value indicates line confusion.}")
    # Footnote entries: literature references and the * (new detection) tag.
    for tgt, letter in sorted(fn_letter.items(), key=lambda kv: kv[1]):
        ref = _latex_escape(lit_refs.get(tgt, {}).get("ref", ""))
        out.append(rf"\tablenotetext{{{letter}}}{{First reported by {ref}.}}")
    if star_used:
        out.append(r"\tablenotetext{*}{First reported in this work.}")
    out.append(r"\end{deluxetable}")
    OUT_TEX.write_text("\n".join(out) + "\n")


def main():
    df = collect()
    print(f"{len(df)} rows assembled")
    write_tex(df)
    print(f"wrote {OUT_TEX}")
    print()
    print(df[["source", "theta_au", "sigma_native", "rrl", "nacl", "kcl", "h2o", "com"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
