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
    "MonR2-IRS3":         "I06059",
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
LIT_MM_NAME = {
    ("NGC6334I", 5):     "NGC6334I-MM1B",
    ("NGC6334I", 4):     "NGC6334I-MM1D",
    ("NGC6334I", 6):     "NGC6334I-MM1C",
    ("NGC6334I", 7):     "NGC6334I-MM1A",
    ("NGC6334I", 2):     "NGC6334I-MM2",
    ("NGC6334I", 1):     "NGC6334I-MM3",
    ("NGC6334IN", 1):    "NGC6334IN-SMA6",
    ("I16547-4247", 1):  "I16547A",
    ("I16547-4247", 2):  "I16547B",
    ("Orion-SrcI", 1):   "Orion-SrcI",
    ("Orion-BN", 1):     "Orion-BN",
    ("OrionBN", 1):      "Orion-BN",
    ("MonR2-IRS3", 1):   "MonR2-IRS3",
    ("MonR2-IRS2", 1):   "MonR2-IRS2-A",
}


def short_handle(name: str, alt_names: str) -> str:
    """A short token usable as the mm-naming prefix.
    IRAS shorthand if known, else the bare target name without the leading
    G... galactic-coord designation, else the name itself."""
    iras = iras_shorthand(name, alt_names)
    if iras and iras != name:
        return iras
    return name


def source_label(target: str, src_id: int, alt_names: str,
                  brightness_rank: int) -> str:
    """Merged Source identifier for the data-summary table.
    Priority: literature mm name → common name + mm{rank} → IRASmm{rank}
    → targetmm{rank}. Common names come from build_target_table.COMMON_NAME.
    Examples: NGC6334I src05 → NGC6334I-MM1B; OrionB-Flame src01 →
    NGC 2024 (Flame)mm1; G268.4222-00.8490 src01 → I09002mm1.
    """
    key = (target, int(src_id))
    if key in LIT_MM_NAME:
        return LIT_MM_NAME[key]
    # Prefer common name if known
    try:
        from build_target_table import COMMON_NAME
    except ImportError:
        COMMON_NAME = {}
    common = COMMON_NAME.get(target)
    if common and common != target:
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
    """Return dict with rrl, salt, com booleans evaluated AT the brightest
    mm continuum source (snr >= 5 at source==bid). Returns None if measurements
    are missing.
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
    salt = bool(sub[sub["group"].isin(["NaCl", "KCl", "H2O"])].shape[0])
    com_pat = re.compile(r"^(CH3OH|CH3CN|HC3N|CH3OCHO|C2H5CN|CH3OCH3|NH2CHO|HCOOH|HNCO)",
                         re.IGNORECASE)
    com = bool(sub[sub["line"].astype(str).apply(lambda L: bool(com_pat.match(L)))].shape[0])
    return dict(rrl=rrl, salt=salt, com=com)


def collect():
    src = pd.read_csv(SRC_CSV)
    lit_path = ROOT / "data/literature_detections.csv"
    lit_targets = set()
    if lit_path.exists():
        try:
            lit_targets = set(pd.read_csv(lit_path)["target"].astype(str))
        except (KeyError, pd.errors.EmptyDataError):
            lit_targets = set()

    # collect per-target best (target, proposal) by analysis_products
    rows = []
    for _, r in src.iterrows():
        name = r["name"]
        iras = field_label(name, r.get("alma_target_names", ""))
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
            det = detections_at_brightest(ad, bid)
            n_det = 0
            if det is not None:
                n_det = sum(int(det[k]) for k in ("rrl", "salt", "com"))
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
                        det = detections_at_brightest(ad, bid)
                        n_det = sum(int(det[k]) for k in ("rrl","salt","com")) if det else 0
                        if (n_det, -theta_au) > (best_key[0], -best_key[1]):
                            best_key = (n_det, theta_au)
                            best = (ad, bid, beam_arcsec, theta_au)
                    break
        if best is None:
            # Still nothing: decide between
            #   (a) target listed in Table 3 (no ALMA <500 AU) → omit entirely
            #   (b) data exist on disk but pipeline not yet run → emit TODO
            n_proposals_500au = 0
            alma_props = str(r.get("alma_proposals", "") or "")
            if alma_props and alma_props != "nan":
                br = r.get("best_res_arcsec")
                if pd.notna(br) and br * dist * 1000.0 < 500.0:
                    n_proposals_500au = 1
            if n_proposals_500au == 0:
                continue
            rows.append(dict(
                source=merged_src_pre,
                theta_au=TODO, sigma_native=TODO, sigma_200au=TODO,
                sigma_200au_10kms=TODO, f5sig=TODO, vwidth=TODO,
                d=f"{dist:.2f}",
                rrl=WIP, salt=WIP, com=WIP,
            ))
            continue
        ad, bid, beam_arcsec, theta_au = best
        rmap = brightness_rank_map(ad / "continuum_sources.csv")
        rank = rmap.get(int(bid), 1)
        merged_src = source_label(name, bid, r.get("alma_target_names", ""), rank)
        sigma_nat, dv_nat = native_noise_for(ad, bid)
        vw = vwidth_for(ad, bid)
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
        rrl = "y" if (det and det["rrl"]) else "n"
        salt = "y" if (det and det["salt"]) else "n"
        com = "y" if (det and det["com"]) else "n"
        if det is None:
            rrl = salt = com = WIP
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
        rows.append(dict(
            source=merged_src,
            theta_au=f"{theta_au:.0f}" if (theta_au and np.isfinite(theta_au)) else TODO,
            sigma_native=f"{sigma_nat*1000.0:.1f}" if sigma_nat else TODO,
            sigma_200au=f"{sigma_200*1000.0:.1f}" if sigma_200 else TODO,
            sigma_200au_10kms=f"{sigma_200_10*1000.0:.2f}" if sigma_200_10 else TODO,
            f5sig=f5_str,
            d=f"{dist:.2f}",
            rrl=rrl, salt=salt, com=com,
        ))
    return pd.DataFrame(rows)


def write_tex(df: pd.DataFrame):
    out = []
    out.append(r"\startlongtable")
    out.append(r"\begin{deluxetable}{lcccccccccc}")
    out.append(r"\tabletypesize{\scriptsize}")
    out.append(r"\tablecaption{Observation summary and detection inventory across the "
               r"full target sample. Beam in AU at the source distance, $\sigma$ in mK. "
               r"Detection columns: y = $\geq 5\sigma$, n = non-detection, "
               r"\textsc{wip}\ = analysis in progress. Source identifier is the "
               r"literature mm/SMA name when known, else "
               r"\texttt{<IRAS>mm<rank>}, ranking by mm peak brightness within "
               r"the field.\label{tab:obs}}")
    out.append(r"\tablehead{")
    out.append(r"\colhead{Source} & "
               r"\colhead{$\theta_\mathrm{maj}$} & "
               r"\colhead{$\sigma_\mathrm{nat}$} & "
               r"\colhead{$\sigma_{200\mathrm{au}}$} & "
               r"\colhead{$\sigma_{200,10}$} & "
               r"\colhead{$f({>}5\sigma)$} & "
               r"\colhead{$d$} & "
               r"\colhead{RRL} & \colhead{NaCl/} & \colhead{KCl} & \colhead{COM} \\")
    out.append(r" & (AU) & (mK) & (mK) & (mK) & & (kpc) & & H$_2$O & & }")
    out.append(r"\startdata")
    for _, r in df.iterrows():
        out.append(
            f"{r['source'].replace('_', r'\\_')} & "
            f"{r['theta_au']} & {r['sigma_native']} & "
            f"{r['sigma_200au']} & {r['sigma_200au_10kms']} & "
            f"{r['f5sig']} & {r['d']} & "
            f"{r['rrl']} & {r['salt']} & n & {r['com']} \\\\"
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
    out.append(r"\end{deluxetable}")
    OUT_TEX.write_text("\n".join(out) + "\n")


def main():
    df = collect()
    print(f"{len(df)} rows assembled")
    write_tex(df)
    print(f"wrote {OUT_TEX}")
    print()
    print(df[["source", "theta_au", "sigma_native", "rrl", "salt"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
