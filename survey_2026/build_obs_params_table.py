"""Fundamental observational parameters table.

For each target with completed line_pipeline analysis:
  - Field (common name + IRAS shorthand)
  - ALMA pointing (RA, Dec from the brightest mm cont source pixel; uncertainty
    estimated from beam_FWHM/(2*SNR))
  - Source velocity (peak velocity of the strongest detected line at the
    brightest source; vlsr_from_literature.json reference if undetected)
  - Reference line used to anchor the velocity (line name with highest SNR
    at the brightest source) — pretty-formatted for LaTeX
  - Detection peak T (mK) for:
      RRL (Hn alpha/beta/gamma)
      H2O 232 GHz
      Brightest NaCl line
    \nodata if not in spectral coverage, "<3 sigma" if covered but undetected.

Outputs:
  data/obs_params.csv
  /orange/adamginsburg/salt/demography_2026/obs_params.tex
"""
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
UVDIR = ROOT / "uvdata"
SRC_CSV = ROOT / "data/sources_L4_d2.csv"
LIT_VLSR = ROOT / "data/vlsr_from_literature.json"
DATA_VLSR = ROOT / "data/vlsr_from_data.json"
OUT_CSV = ROOT / "data/obs_params.csv"
OUT_TEX = Path("/orange/adamginsburg/salt/demography_2026/obs_params.tex")

IRAS_RE = re.compile(r"(?:IRAS[ _]*|\bI)([0-9]{5})(?:[+-][0-9]+)?", re.IGNORECASE)


def field_label(name, alt):
    """Compact field display name for Table 5: prefer COMMON_NAME (so the
    label matches Tables 1/4), else IRAS shorthand from alma_target_names
    (so T5 matches T1's `IRAS XXXXX+YYYY` fallback), else the bare target."""
    try:
        from build_target_table import COMMON_NAME, alma_target_names_to_iras
    except ImportError:
        COMMON_NAME = {}
        alma_target_names_to_iras = lambda _s: None
    common = COMMON_NAME.get(name)
    if common:
        return common
    iras = alma_target_names_to_iras(alt)
    if iras:
        return iras
    return name


def brightest_source_row(proposal_dir):
    cont = proposal_dir / "continuum_sources.csv"
    if not cont.exists() or cont.stat().st_size == 0:
        return None
    try:
        df = pd.read_csv(cont)
    except pd.errors.EmptyDataError:
        return None
    if df.empty or "peak_Jybeam" not in df.columns:
        return None
    idx = df["peak_Jybeam"].idxmax()
    return df.loc[idx]


def beam_arcsec_from_first_cube(proposal_dir):
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
        if h.get("BMAJ") is not None:
            return float(h["BMAJ"]) * 3600.0
    return None


def best_proposal(name, dist):
    """Return (proposal_dir, brightest_id, beam_arcsec) for the proposal with
    the most lines detected at its own brightest source (tiebreak: smallest
    beam). Mirrors build_data_summary."""
    adirs = sorted((ANALYSIS / name).glob("2*")) if (ANALYSIS / name).is_dir() else []
    if not adirs:
        return None
    best = None
    best_key = (-1, np.inf)
    for ad in adirs:
        row = brightest_source_row(ad)
        if row is None:
            continue
        bid = int(row["id"])
        meas = ad / "line_measurements.csv"
        n_det = 0
        if meas.exists():
            try:
                m = pd.read_csv(meas)
                if not m.empty and "snr" in m.columns:
                    n_det = int(((m["source"] == bid) & (m["snr"] >= 5.0)).sum())
            except pd.errors.EmptyDataError:
                pass
        beam = beam_arcsec_from_first_cube(ad)
        theta_au = beam * dist * 1000.0 if beam else np.inf
        key = (n_det, -theta_au)
        if key > (best_key[0], -best_key[1]):
            best_key = (n_det, theta_au)
            best = (ad, bid, beam, row)
    return best


RRL_RE = re.compile(r"^H\d+(alpha|beta|gamma|delta)$")
H2O_RE = re.compile(r"^H2O.*232")
NACL_RE = re.compile(r"^NaCl_")
KCL_RE = re.compile(r"^KCl_")


def load_lit_refs():
    """Match build_data_summary.load_lit_refs (kept inline to avoid an import
    cycle); returns {target: {ref, species_kind}} including alias spellings."""
    lit_path = ROOT / "data/literature_detections.csv"
    if not lit_path.exists():
        return {}
    try:
        df = pd.read_csv(lit_path)
    except pd.errors.EmptyDataError:
        return {}
    species = ["NaCl", "KCl", "H2O", "SiO", "SiS", "SO", "COM", "RRL"]
    out = {}
    for _, r in df.iterrows():
        target = str(r["target"])
        kinds = {s: str(r.get(f"{s}_kind", "") or "").strip().lower()
                  for s in species}
        if any(v in ("det", "tent") for v in kinds.values()):
            entry = dict(ref=str(r["reference"]), species_kind=kinds)
            out[target] = entry
            for alias in (target.replace("-", "_"), target.replace("_", "-")):
                if alias != target:
                    out.setdefault(alias, entry)
    return out


def deg_to_alma_str(ra_deg, dec_deg):
    """ALMA-style colon-separated J2000 string."""
    import astropy.coordinates as ac
    import astropy.units as u
    c = ac.SkyCoord(ra_deg * u.deg, dec_deg * u.deg)
    ra = c.ra.to_string(unit=u.hourangle, sep=":", pad=True, precision=2)
    dec = c.dec.to_string(unit=u.deg, sep=":", pad=True, precision=1, alwayssign=True)
    return ra, dec


def pos_uncert_arcsec(beam_arcsec, snr):
    """Position uncertainty ~ beam_FWHM / (2 * SNR), Reid+1988 / synth-array
    rule of thumb."""
    if not beam_arcsec or not snr or snr <= 0:
        return None
    return float(beam_arcsec) / (2.0 * float(snr))


def fmt_T(val_K):
    if val_K is None or not np.isfinite(val_K):
        return r"\nodata"
    if abs(val_K) >= 1.0:
        return f"{val_K:.2f}"
    return f"{val_K*1000:.0f}"  # in mK


def fmt_kkms(val):
    if val is None or not np.isfinite(val):
        return r"\nodata"
    if abs(val) >= 1.0:
        return f"{val:.2f}"
    return f"{val*1000:.0f}"  # in mK km/s


def line_value(meas_df, bid, regex, integ_unit):
    """Return (peak_T_K, integ_K_kms, line_name, snr, status, sigma_K).
    sigma_K is the per-channel noise (K) used to compute 3σ upper limits
    in cells with status='covered'.
    """
    matches = meas_df[meas_df["line"].astype(str).apply(lambda L: bool(regex.match(L)))]
    any_in_band = matches.shape[0] > 0
    if not any_in_band:
        return None, None, None, None, "no", None
    at_src = matches[matches["source"] == bid]
    if at_src.empty:
        # in band but not measured at brightest source; report 3σ at the
        # tightest noise across any source covered by the same line set.
        sig = float(np.nanmin(matches["sigma"])) if "sigma" in matches.columns else None
        return None, None, None, None, "covered", sig
    best = at_src.nlargest(1, "snr").iloc[0]
    snr = float(best["snr"])
    sig = float(best["sigma"]) if "sigma" in at_src.columns else None
    if snr < 5.0:
        return None, None, None, snr, "covered", sig
    return (
        float(best["peak_Kkms_or_unit"]),  # ambiguous name; treated as peak T
        float(best["integ"]),
        str(best["line"]),
        snr,
        "detected",
        sig,
    )


_LIT_VLSR = json.loads(LIT_VLSR.read_text()) if LIT_VLSR.exists() else {}
_DATA_VLSR = json.loads(DATA_VLSR.read_text()) if DATA_VLSR.exists() else {}


def latex_line(line_name):
    """Format an internal line label into compact LaTeX. Examples:
       H30alpha       → H30$\\alpha$
       H2O_v2_5_5_...232 → H$_2$O ($v_2$, 232\\,GHz)
       NaCl_v0_J18-17 → NaCl $J$=18--17 ($v$=0)
       SiO_v0_5-4     → SiO $J$=5--4
       CH3OH_4-3      → CH$_3$OH $J$=4--3
       SO_65-54       → SO
    """
    if not isinstance(line_name, str) or not line_name:
        return r"\nodata"
    s = line_name
    # RRL
    m = re.match(r"^H(\d+)(alpha|beta|gamma|delta)$", s)
    if m:
        greek = {"alpha": r"\alpha", "beta": r"\beta",
                  "gamma": r"\gamma", "delta": r"\delta"}[m.group(2)]
        return f"H{m.group(1)}${greek}$"
    # H2O 232 GHz (any quantum-number suffix)
    m = re.match(r"^H2O(?:_v(\d))?", s)
    if m:
        if m.group(1):
            return f"H$_2$O ($v_{{{m.group(1)}}}$, 232\\,GHz)"
        return r"H$_2$O 232\,GHz"
    # NaCl / KCl with optional v=N and J=A-B
    m = re.match(r"^(NaCl|KCl|Na37Cl|K37Cl)(?:_v(\d))?(?:_J(\d+)-(\d+))?", s)
    if m:
        mol = m.group(1).replace("37", r"$^{37}$")
        parts = [mol]
        if m.group(3):
            parts.append(rf" $J$={m.group(3)}--{m.group(4)}")
        if m.group(2):
            parts.append(rf" ($v$={m.group(2)})")
        return "".join(parts)
    # SiO/SiS with v and J
    m = re.match(r"^(SiO|SiS)(?:_v(\d))?(?:_J?(\d+)-(\d+))?", s)
    if m:
        out = m.group(1)
        if m.group(3):
            out += rf" $J$={m.group(3)}--{m.group(4)}"
        return out
    # SO / SO2
    if s.startswith("SO2"):
        return r"SO$_2$"
    if s.startswith("SO"):
        m = re.match(r"^SO_(\d+)_(\d+)-(\d+)_(\d+)", s)
        if m:
            return rf"SO $J_{{N}}$={m.group(1)}$_{{{m.group(2)}}}$--{m.group(3)}$_{{{m.group(4)}}}$"
        return "SO"
    # COMs (simple molecule label, drop transition)
    for mol_raw, mol_tex in [("CH3OCHO", r"CH$_3$OCHO"), ("CH3OCH3", r"CH$_3$OCH$_3$"),
                              ("CH3OH", r"CH$_3$OH"), ("CH3CN", r"CH$_3$CN"),
                              ("C2H5CN", r"C$_2$H$_5$CN"), ("HC3N", r"HC$_3$N"),
                              ("NH2CHO", r"NH$_2$CHO"), ("HCOOH", r"HCOOH"),
                              ("HNCO", r"HNCO"), ("13CO", r"$^{13}$CO"),
                              ("C18O", r"C$^{18}$O"), ("H2CO", r"H$_2$CO")]:
        if s.startswith(mol_raw):
            return mol_tex
    # Fallback: escape underscores
    return s.replace("_", r"\_")


def collect():
    src = pd.read_csv(SRC_CSV)
    rows = []
    for _, r in src.iterrows():
        name = r["name"]
        dist = float(r["dist_kpc"])
        bp = best_proposal(name, dist)
        if bp is None:
            continue
        ad, bid, beam_arcsec, srow = bp
        meas_csv = ad / "line_measurements.csv"
        if not meas_csv.exists():
            continue
        try:
            meas = pd.read_csv(meas_csv)
        except pd.errors.EmptyDataError:
            continue
        if meas.empty or not {"source", "snr", "line", "peak_v"} <= set(meas.columns):
            continue
        ra_str, dec_str = deg_to_alma_str(float(srow["ra_deg"]), float(srow["dec_deg"]))
        snr_cont = float(srow.get("snr", 0.0)) if "snr" in srow.index else 0.0
        sigpos = pos_uncert_arcsec(beam_arcsec, snr_cont)
        # source velocity: pick highest-SNR detected line at brightest src;
        # use its peak_v. If nothing detected, fall back to the
        # vlsr_from_literature reference (no "canonical" launch value).
        at_src = meas[(meas["source"] == bid) & (meas["snr"] >= 5.0)]
        if not at_src.empty:
            best_line = at_src.nlargest(1, "snr").iloc[0]
            vsrc = float(best_line["peak_v"])
            vref_line = str(best_line["line"])
        else:
            vlit = _LIT_VLSR.get(name) or _DATA_VLSR.get(name)
            if vlit and vlit.get("v_LSR_kms") is not None:
                vsrc = float(vlit["v_LSR_kms"])
                ref_short = str(vlit.get("reference", "literature")).split(";")[0][:32]
                vref_line = f"lit: {ref_short}"
            else:
                vsrc = None
                vref_line = ""
        # line values
        rrl_peak, rrl_int, rrl_line, _, rrl_status, rrl_sig = line_value(meas, bid, RRL_RE, "")
        h2o_peak, h2o_int, h2o_line, _, h2o_status, h2o_sig = line_value(meas, bid, H2O_RE, "")
        nacl_peak, nacl_int, nacl_line, _, nacl_status, nacl_sig = line_value(meas, bid, NACL_RE, "")
        kcl_peak, kcl_int, kcl_line, _, kcl_status, kcl_sig = line_value(meas, bid, KCL_RE, "")
        rows.append(dict(
            name=name,
            field=field_label(name, r.get("alma_target_names", "")),
            ra=ra_str, dec=dec_str,
            pos_unc_arcsec=sigpos,
            vsrc_kms=vsrc, vref_line=vref_line,
            rrl_line=rrl_line, rrl_peak_K=rrl_peak, rrl_int=rrl_int, rrl_status=rrl_status, rrl_sig_K=rrl_sig,
            h2o_line=h2o_line, h2o_peak_K=h2o_peak, h2o_int=h2o_int, h2o_status=h2o_status, h2o_sig_K=h2o_sig,
            nacl_line=nacl_line, nacl_peak_K=nacl_peak, nacl_int=nacl_int, nacl_status=nacl_status, nacl_sig_K=nacl_sig,
            kcl_line=kcl_line, kcl_peak_K=kcl_peak, kcl_int=kcl_int, kcl_status=kcl_status, kcl_sig_K=kcl_sig,
            beam_arcsec=beam_arcsec, source_id=bid,
        ))
    # Dedupe by (ra, dec, field) within 1": two source-CSV rows that resolve
    # to the same continuum source (e.g. G336.4917A and G336.4917B both
    # → I16362) collapse to one table row, keeping the first.
    if rows:
        deduped = []
        seen = []
        for row in rows:
            key_ra, key_dec = row["ra"], row["dec"]
            if any(s == (key_ra, key_dec) for s in seen):
                continue
            seen.append((key_ra, key_dec))
            deduped.append(row)
        rows = deduped
    return pd.DataFrame(rows)


def fmt_cell(peak, integ, status, lit_kind=None, lit_letter=None, sigma_K=None):
    """One cell. Detection: peak in mK. Covered non-detection: explicit 3σ
    upper limit `<N` in mK. No coverage: \\nodata. Literature-detected
    cells with no re-measurement are decorated with the lit footnote.
    """
    if status == "detected":
        if peak is None or not np.isfinite(peak):
            return r"\nodata"
        return f"{peak*1000:.1f}"
    if status == "covered":
        if sigma_K is not None and np.isfinite(sigma_K):
            ul_mK = 3.0 * sigma_K * 1000.0
            base = f"$<{ul_mK:.1f}$"
        else:
            base = r"$<\!3\sigma$"
    else:
        base = r"\nodata"
    if lit_kind in ("det", "tent") and lit_letter:
        return rf"{base}\,\tablenotemark{{{lit_letter}}}"
    return base


def lit_ref_to_num(text, citet_map, footnotes):
    """Replace 'lit: <free text>' with 'lit$^{[N]}$' and register the
    free-text reference in footnotes/citet_map."""
    if not text or not text.startswith("lit:"):
        return text
    raw = text[4:].strip()
    if not raw:
        return "lit"
    # Escape LaTeX-special chars before storing in footnote
    safe = raw
    for ch, esc in [("&", r"\&"), ("_", r"\_"), ("%", r"\%"),
                     ("#", r"\#"), ("~", r"$\sim$")]:
        safe = safe.replace(ch, esc)
    if safe in citet_map:
        n = citet_map[safe]
    else:
        n = len(citet_map) + 1
        citet_map[safe] = n
        footnotes[n] = safe
    return rf"lit$^{{[{n}]}}$"


def write_tex(df):
    df = df.sort_values("name")
    lit_refs = load_lit_refs()
    citet_map = {}
    footnotes = {}
    # Literature footnote letters per target (for cells)
    fn_letter = {}
    out = []
    out.append(r"\begin{deluxetable}{lcccccccccc}")
    out.append(r"\tabletypesize{\tiny}")
    out.append(r"\tablecaption{Fundamental observational parameters per target. "
               r"Coordinates are the J2000 position of the brightest mm "
               r"continuum source within the analyzed ALMA field; the "
               r"positional uncertainty is the synthesized-beam major axis "
               r"divided by twice the source signal-to-noise ratio. "
               r"$v_\mathrm{LSR}$ is the peak velocity of the highest-SNR line "
               r"detected at the brightest source (Ref.\ line column); for "
               r"non-detections the literature value is reported (Ref.\ "
               r"line = `lit: \dots'). Detection cells give peak brightness "
               r"(mK) of the highest-SNR line in the species class. "
               r"`$<\!N$' in a covered cell is the $3\sigma$ upper limit "
               r"in mK at the brightest source's peak pixel; "
               r"\nodata\ marks no spectral coverage of the species. "
               r"Species-cell footnote markers cite literature detections "
               r"for which we do not have a re-measured value here.\label{tab:obspars}}")
    out.append(r"\tablehead{")
    out.append(r"\colhead{Field} & \colhead{R.A.} & \colhead{Dec.} & "
               r"\colhead{$\sigma_\mathrm{pos}$} & "
               r"\colhead{$v_\mathrm{LSR}$} & \colhead{Ref.\ line} & "
               r"\colhead{RRL} & \colhead{H$_2$O\,232} & "
               r"\colhead{NaCl} & \colhead{KCl} \\")
    out.append(r" & (J2000) & (J2000) & ($''$) & (\kms) & & "
               r"(mK) & (mK) & (mK) & (mK) }")
    out.append(r"\startdata")
    for _, r in df.iterrows():
        field = str(r["field"]).replace("_", r"\_")
        ref_raw = str(r["vref_line"])
        if ref_raw.startswith("lit:"):
            ref = lit_ref_to_num(ref_raw, citet_map, footnotes)
        elif ref_raw:
            ref = latex_line(ref_raw)
        else:
            ref = r"\nodata"
        # Literature footnote letter for this target (if any lit det/tent)
        tgt = r["name"]
        lit_entry = lit_refs.get(tgt) or lit_refs.get(tgt.replace("-", "_")) \
                    or lit_refs.get(tgt.replace("_", "-"))
        lit_kinds = lit_entry["species_kind"] if lit_entry else {}
        # Allocate a footnote letter only when at least one species cell will
        # actually consume it (lit says det/tent AND we have no detection).
        needs_mark = False
        for sp, status_col in (("NaCl", "nacl_status"), ("KCl", "kcl_status"),
                                ("H2O", "h2o_status"), ("RRL", "rrl_status")):
            if lit_kinds.get(sp) in ("det", "tent") and r[status_col] != "detected":
                needs_mark = True
                break
        lit_letter = None
        if needs_mark:
            if tgt not in fn_letter:
                fn_letter[tgt] = chr(ord('a') + len(fn_letter))
            lit_letter = fn_letter[tgt]
        rrl_l = lit_kinds.get("RRL")
        h2o_l = lit_kinds.get("H2O")
        nacl_l = lit_kinds.get("NaCl")
        kcl_l = lit_kinds.get("KCl")
        rrl = fmt_cell(r["rrl_peak_K"], r["rrl_int"], r["rrl_status"], rrl_l, lit_letter, r["rrl_sig_K"])
        h2o = fmt_cell(r["h2o_peak_K"], r["h2o_int"], r["h2o_status"], h2o_l, lit_letter, r["h2o_sig_K"])
        nacl = fmt_cell(r["nacl_peak_K"], r["nacl_int"], r["nacl_status"], nacl_l, lit_letter, r["nacl_sig_K"])
        kcl = fmt_cell(r["kcl_peak_K"], r["kcl_int"], r["kcl_status"], kcl_l, lit_letter, r["kcl_sig_K"])
        sig_str = f"{r['pos_unc_arcsec']:.3f}" if pd.notna(r["pos_unc_arcsec"]) else r"\nodata"
        v_str = f"{r['vsrc_kms']:+.1f}" if pd.notna(r["vsrc_kms"]) else r"\nodata"
        out.append(
            f"{field} & {r['ra']} & {r['dec']} & {sig_str} & "
            f"{v_str} & {ref} & {rrl} & {h2o} & {nacl} & {kcl} \\\\"
        )
    out.append(r"\enddata")
    if footnotes:
        ref_strs = [rf"[{n}] {footnotes[n]}" for n in sorted(footnotes)]
        out.append(r"\tablerefs{" + "; ".join(ref_strs) + "}")
    for tgt, letter in sorted(fn_letter.items(), key=lambda kv: kv[1]):
        ref = lit_refs[tgt]["ref"].replace("&", r"\&").replace("_", r"\_")
        out.append(rf"\tablenotetext{{{letter}}}{{Literature detection: {ref}.}}")
    out.append(r"\end{deluxetable}")
    OUT_TEX.write_text("\n".join(out) + "\n")


def main():
    df = collect()
    print(f"{len(df)} rows")
    if df.empty:
        return
    df.to_csv(OUT_CSV, index=False)
    write_tex(df)
    print(f"wrote {OUT_TEX}")
    print(df[["name", "ra", "dec", "pos_unc_arcsec", "vsrc_kms",
              "rrl_status", "h2o_status", "nacl_status"]].to_string(index=False))


if __name__ == "__main__":
    main()
