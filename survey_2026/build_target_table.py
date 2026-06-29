"""Build LaTeX target table for the demography_2026 paper.

For each source in sources_L4_d2.csv, query ALMA archive for project codes
whose spatial_resolution_arcsec * dist_pc < 500 AU; emit a longtable to
/orange/adamginsburg/salt/demography_2026/targets.tex.

Columns: Name, RA, Dec, d [kpc], Lbol [Lsun], ALMA proposals (<500 AU).
Distance + luminosity reference citation keys assigned per origin catalog.
"""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import astropy.units as u
import astropy.coordinates as ac
from astroquery.alma import Alma

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
SRC = ROOT / "data/sources_L4_d2.csv"
OUT = Path("/orange/adamginsburg/salt/demography_2026/targets.tex")

# Refs by origin catalog. Special cases override per name.
REFS_BY_ORIGIN = {
    "RMS":   r"\citet{Lumsden2013,Urquhart2014}",
    "HiGAL": r"\citet{Elia2017}",
}
LUM_REF_BY_ORIGIN = {
    "RMS":   r"\citet{Mottram2011,Lumsden2013}",
    "HiGAL": r"\citet{Elia2017}",
}
SPECIAL_DIST_REF = {
    "OrionSrcI":           r"\citet{Menten2007}",
    "Orion_SrcI":          r"\citet{Menten2007}",
    "Orion-SrcI":          r"\citet{Menten2007}",
    "Orion-BN":            r"\citet{Menten2007}",
    "OrionBN":             r"\citet{Menten2007}",
    "OrionB-Flame":        r"\citet{Kounkel2017}",
    "S140-IRS1":           r"\citet{Hirota2008}",
    "NGC6334I":            r"\citet{Chibueze2014}",
    "NGC6334IN":           r"\citet{Chibueze2014}",
    "IRAS17233-3606":      r"\citet{Chibueze2014}",
    "MonR2-IRS3":          r"\citet{Herbst1976}",
    "MonR2-IRS2":          r"\citet{Herbst1976}",
    "GGD12-15":            r"\citet{Gomez2002}",
    "NGC6514":             r"\citet{Kuhn2019}",
    "Lagoon-Her36":        r"\citet{Damiani2019}",
    "G353.2+0.9":          r"\citet{Kuhn2019}",
    "G268.4222-00.8490":   r"\citet{Getman2019}",
    "G010.8411-02.5919":   r"\citet{Zucker2020}",
    "G081.6802+00.5405B":  r"\citet{Rygl2012}",
    "G081.6802+00.5405A":  r"\citet{Rygl2012}",
    "G345.0052+01.8209":   r"\citet{Lopez2011}",
    "G345.0061+01.7944C":  r"\citet{Lopez2011}",
    "G345.2244+01.0304":   r"\citet{Lopez2011}",
    "G345.4938+01.4677":   r"\citet{Guzman2014}",
    "G345.5043+00.3480":   r"\citet{Urquhart2018}",
}
SPECIAL_LUM_REF = {
    "OrionSrcI":          r"\citet{Testi2010}",
    "Orion_SrcI":         r"\citet{Testi2010}",
    "Orion-SrcI":         r"\citet{Testi2010}",
    "Orion-BN":           r"\citet{Gezari1998}",
    "OrionBN":            r"\citet{Gezari1998}",
    "OrionB-Flame":       r"\citet{Bik2003}",
    "S140-IRS1":          r"\citet{Maud2013}",
    "NGC6334I":           r"\citet{Sandell2000}",
    "NGC6334IN":          r"\citet{Sandell2000}",
    "IRAS17233-3606":     r"\citet{Leurini2008}",
    "MonR2-IRS3":         r"\citet{Henning1992}",
    "MonR2-IRS2":         r"\citet{Henning1992}",
    "GGD12-15":           r"\citet{Gomez2002}",
    "NGC6514":            r"\citet{Rho2006}",
    "Lagoon-Her36":       r"\citet{Goto2006,Arias2010}",
    "G353.2+0.9":         r"\citet{Giannetti2012}",
    "G268.4222-00.8490":  r"\citet{Getman2019}",
}


def query_alma_lt500au(coord, dist_kpc, radius=15 * u.arcsec):
    """Return (per_code, best_res_au, n_obs_total).

    per_code: list of (code, min_res_arcsec) sorted by code, restricted to obs
              whose spatial resolution is finer than 500 AU at the source d.
    best_res_au: best resolution (smallest) at this position over ALL obs (AU).
    n_obs_total: total number of obs returned by ALMA archive at this position.
    """
    try:
        q = Alma.query_region(coord, radius=radius)
    except (ConnectionError, ValueError) as e:
        print(f"  query error: {e}")
        return [], None, 0
    if q is None or len(q) == 0:
        return [], None, 0
    sr_col = "spatial_resolution" if "spatial_resolution" in q.colnames else "s_resolution"
    pc_col = "proposal_id" if "proposal_id" in q.colnames else "project_code"
    sr = np.asarray(q[sr_col], dtype=float)
    pc = np.asarray(q[pc_col], dtype=str)
    finite = np.isfinite(sr) & (sr > 0)
    n_obs_total = int(finite.sum())
    if not finite.any():
        return [], None, 0
    best_res_arcsec = float(np.nanmin(sr[finite]))
    best_res_au = best_res_arcsec * dist_kpc * 1000.0
    res_au = sr * dist_kpc * 1000.0
    mask = finite & (res_au < 500.0)
    # Group by code → min resolution
    per_code = {}
    for c, s in zip(pc[mask], sr[mask]):
        per_code[c] = min(per_code.get(c, np.inf), float(s))
    return sorted(per_code.items()), best_res_au, n_obs_total


def fmt_coord(ra_deg, dec_deg):
    """Return a single merged J2000 coord string `HH:MM:SS.s±DD:MM:SS` with
    colons but no whitespace between RA and Dec (Dec sign separates them)."""
    c = ac.SkyCoord(ra_deg * u.deg, dec_deg * u.deg)
    ra = c.ra.to_string(unit=u.hourangle, sep=":", pad=True, precision=1)
    dec = c.dec.to_string(unit=u.deg, sep=":", pad=True, precision=0, alwayssign=True)
    return f"{ra}{dec}"


def fmt_galactic(glon, glat):
    """Compact Gxxx.xxxx{+|-}xx.xxxx wrapped in \\texttt{} so the +/- glyphs
    align in a fixed-width font."""
    s = "+" if glat >= 0 else "-"
    return rf"\texttt{{G{float(glon):08.4f}{s}{abs(float(glat)):07.4f}}}"


def fmt_distance(d_kpc):
    """Distance formatted to appropriate sig figs."""
    if d_kpc >= 10:
        return f"{d_kpc:.0f}"
    if d_kpc >= 1:
        return f"{d_kpc:.1f}"
    return f"{d_kpc:.2f}"


# Common / IRAS / popular name overrides for table display.
COMMON_NAME = {
    "OrionSrcI":          "Orion-SrcI",
    "Orion_SrcI":         "Orion-SrcI",
    "OrionBN":            "Orion-BN",
    "Orion-BN":           "Orion-BN",
    "OrionBN-KL":         "Orion-BN/KL",
    "OrionB-Flame":       "NGC 2024 (Flame)",
    "NGC6334I":           "NGC 6334I-mm1b",
    "NGC6334IN":          "NGC 6334I(N)-SMA6",
    "MonR2-IRS3":         "MonR2-IRS3",
    "MonR2-IRS2":         "MonR2-IRS2",
    "GGD12-15":           "GGD 12-15",
    "Lagoon-Her36":       "Her 36 (M8)",
    "NGC6514":            "NGC 6514 (M20)",
    "S140-IRS1":          "S140-IRS1",
    "G268.4222-00.8490":  "G268.42-0.85",
    "G010.8411-02.5919":  "GGD 27 / I18162",
    "G011.9197-00.6131":  "G11.92 mm1",
    "G017.6396+00.1580":  "G17.64+0.16",
    "G351.7745-00.5377":  "G351.77 mm1",
    "G013.6562-00.5997":  "W33A",
    "G019.6097-00.2342":  "I18243-1210",
    "G023.0099-00.4108":  "G23.01 mm2",
    "G028.8621+00.0656":  "G28.86 mm1",
    "G029.0036+00.0083":  "G29.00 mm1",
    "G033.6437-00.2277":  "K3-50A",
    "G035.5781+00.0709":  "I18553+0214",
    "G045.0712+00.1321":  "K47",
    "G081.6802+00.5405A": "DR21(OH) A",
    "G081.6802+00.5405B": "DR21(OH) B",
    "G133.7184+01.2237":  "W3-IRS5",
    "G133.7150+01.2155":  "W3-IRS4",
    "G133.7354+01.1827":  "W3(OH)",
    "G189.7791-00.1037":  "S252A",
    "G192.6005-00.0479":  "S255IR SMA1",
    "G232.6207+00.9959":  "I07299-1651",
    "G268.3957-00.8541":  "G268.40-0.85",
    "G343.5036-00.0119":  "I16957-4306",
    "G345.0061+01.7944C": "G345.01 C",
    "G345.0052+01.8209":  "G345.00+1.82",
    "G345.4938+01.4677":  "IRAS 16562-3959",
    "G345.5043+00.3480":  "I17016-4124",
    "I16547-4247":        "I16547-4247",
}


def alma_target_names_to_iras(s):
    if not isinstance(s, str) or not s.strip():
        return None
    import re as _re
    m = _re.search(r"IRAS\s*([0-9]{5}[+-][0-9]{4})", s)
    if m:
        return f"IRAS {m.group(1)}"
    m = _re.search(r"\bI([0-9]{5})\b", s)
    if m:
        return f"I{m.group(1)}"
    return None


def common_name_for(row):
    if row.name in COMMON_NAME:
        return COMMON_NAME[row.name]
    iras = alma_target_names_to_iras(getattr(row, "alma_target_names", ""))
    if iras:
        return iras
    return row.name


def load_vlsr_sources():
    """Return dict {target: dict(v_kms, method, ref, tracers)}.

    Preference: vlsr_from_data.json > vlsr_from_literature.json > RMS fits.
    """
    out = {}
    lit_path = ROOT / "data/vlsr_from_literature.json"
    if lit_path.exists():
        lit = json.loads(lit_path.read_text())
        for k, v in lit.items():
            if v.get("v_LSR_kms") is None:
                continue
            out[k] = dict(v_kms=float(v["v_LSR_kms"]),
                           method="literature",
                           ref=v.get("reference", ""),
                           tracers=None)
    data_path = ROOT / "data/vlsr_from_data.json"
    if data_path.exists():
        dat = json.loads(data_path.read_text())
        for k, v in dat.items():
            if v.get("v_LSR_kms") is None:
                continue
            tracers = list((v.get("tracers") or {}).keys())
            out[k] = dict(v_kms=float(v["v_LSR_kms"]),
                           method="data",
                           ref="",
                           tracers=tracers)
    # RMS fallback
    rms_path = ROOT / "data/rms.fits"
    if rms_path.exists():
        from astropy.table import Table
        t = Table.read(rms_path)
        names = np.array([str(n) if not isinstance(n, bytes) else n.decode()
                          for n in t["Name"]])
        for k in names:
            if k in out:
                continue
            idx = np.where(names == k)[0]
            if not len(idx):
                continue
            v = float(t["vLSR"][idx[0]])
            if not np.isfinite(v):
                continue
            out[k] = dict(v_kms=v, method="rms",
                           ref=r"\citet{Lumsden2013}", tracers=None)
    return out


def main():
    src = pd.read_csv(SRC).sort_values("dist_kpc")
    vlsr_lookup = load_vlsr_sources()
    rows = []
    footnotes = {}  # marker -> footnote text

    # Numbered \citet footnote ledger: any \citet{...} inside a table cell
    # gets replaced with a [N] superscript number; the actual \citet command
    # goes into a \tablenotetext at the bottom of the table.
    citet_map = {}  # citet text -> numeric tag

    def citet_to_num(text):
        """Replace every \citet{...} in `text` with a [N] superscript;
        register the citet command in citet_map / footnotes."""
        import re as _re
        def repl(m):
            cmd = m.group(0)  # \citet{...}
            if cmd not in citet_map:
                n = len(citet_map) + 1
                citet_map[cmd] = n
            return rf"$^{{[{citet_map[cmd]}]}}$"
        return _re.sub(r"\\citet\{[^}]+\}", repl, text)
    for i, r in enumerate(src.itertuples(), 1):
        coord = ac.SkyCoord(r.ra_deg * u.deg, r.dec_deg * u.deg)
        print(f"[{i}/{len(src)}] {r.name} d={r.dist_kpc:.2f} kpc ...", flush=True)
        per_code, best_res_au, n_obs_total = query_alma_lt500au(coord, r.dist_kpc)
        codes = [c for c, _ in per_code]
        coord = fmt_coord(r.ra_deg, r.dec_deg)
        common = common_name_for(r)
        gname = fmt_galactic(r.glon, r.glat)
        d_str = fmt_distance(r.dist_kpc)
        dref = SPECIAL_DIST_REF.get(r.name, REFS_BY_ORIGIN.get(r.origin, ""))
        lref = SPECIAL_LUM_REF.get(r.name, LUM_REF_BY_ORIGIN.get(r.origin, ""))
        # v_LSR + reference. Long literature/data references go into per-source
        # footnotes; the table cell just carries the method label + a small
        # tablenotemark letter.
        vrec = vlsr_lookup.get(r.name)
        if vrec is None:
            v_str = r"\nodata"
            vref = ""
        else:
            v = vrec["v_kms"]
            method = vrec["method"]
            tag = chr(ord('a') + len(footnotes))
            if method == "data":
                tracers = vrec.get("tracers") or []
                tracers_str = (
                    ", ".join(t.replace("_", r"\_") for t in tracers)
                    if tracers else "(none recorded)"
                )
                safe_name = r.name.replace("_", r"\_")
                footnotes[tag] = (
                    f"Source {safe_name}: $v_\\mathrm{{LSR}}$ measured from "
                    f"in-band corroborating lines: {tracers_str}."
                )
                v_str = f"{v:+.1f}"
                vref = rf"data\tablenotemark{{{tag}}}"
            elif method == "literature":
                safe_name = r.name.replace("_", r"\_")
                ref_txt = vrec.get("ref", "")
                # Escape special LaTeX chars in free-text references
                for ch_from, ch_to in [("&", r"\&"), ("_", r"\_"),
                                        ("%", r"\%"), ("#", r"\#"),
                                        ("~", r"$\sim$")]:
                    ref_txt = ref_txt.replace(ch_from, ch_to)
                if len(ref_txt) > 220:
                    ref_txt = ref_txt[:215] + "..."
                footnotes[tag] = (
                    f"Source {safe_name}: $v_\\mathrm{{LSR}}$ literature "
                    f"reference: {ref_txt}."
                )
                v_str = f"{v:+.1f}"
                vref = rf"lit.\tablenotemark{{{tag}}}"
            else:
                v_str = f"{v:+.1f}"
                vref = vrec.get("ref", "RMS")
        # Convert any \citet{} in references to numbered footnotes
        dref_n = citet_to_num(dref)
        lref_n = citet_to_num(lref)
        vref_n = citet_to_num(vref)
        rows.append(dict(
            name=r.name.replace("_", r"\_"),
            common=common.replace("_", r"\_"),
            gname=gname,
            coord=coord,
            d=r.dist_kpc, d_str=d_str, lbol=r.lbol_lsun,
            dref=dref_n, lref=lref_n,
            v_lsr=v_str, vref=vref_n,
            codes=codes,
            per_code=per_code,
            best_res_au=best_res_au,
            n_obs_total=n_obs_total,
        ))

    pd.DataFrame([{**r, "codes": "; ".join(r["codes"])} for r in rows]).to_csv(
        ROOT / "data/target_table_alma.csv", index=False
    )

    lines = []
    # Compact portrait version of the targets table. Common name + Galactic
    # coord drive identification; J2000 RA/Dec are space-stripped to fit.
    lines.append(r"\startlongtable")
    lines.append(r"\begin{deluxetable}{lllcccccc}")
    lines.append(r"\tabletypesize{\scriptsize}")
    lines.append(r"\tablecaption{Target source list. ALMA project codes are listed when "
                 r"the corresponding observations achieved a spatial resolution finer than "
                 r"500\,AU at the source distance.\label{tab:targets}}")
    lines.append(r"\tablehead{")
    lines.append(r"\colhead{Source} & \colhead{Galactic} & "
                 r"\colhead{R.A.\,Dec.} & "
                 r"\colhead{$d$} & \colhead{$d^{*}$} & "
                 r"\colhead{$L_\mathrm{bol}$} & \colhead{$L^{*}$} & "
                 r"\colhead{$v_\mathrm{LSR}$} & \colhead{$v^{*}$} \\")
    lines.append(r"& (l\,b) & (J2000) & (kpc) & ref. & "
                 r"($10^4\,L_\odot$) & ref. & (km\,s$^{-1}$) & ref.}")
    lines.append(r"\startdata")
    for r in rows:
        lines.append(
            f"{r['common']} & {r['gname']} & \\texttt{{{r['coord']}}} & "
            f"{r['d_str']} & {r['dref']} & "
            f"{r['lbol']/1e4:.2f} & {r['lref']} & "
            f"{r['v_lsr']} & {r['vref']} \\\\"
        )
    lines.append(r"\enddata")
    # Numbered citet references
    refs_list = sorted(citet_map.items(), key=lambda kv: kv[1])
    if refs_list:
        ref_strs = [rf"[{n}] {cmd}" for cmd, n in refs_list]
        lines.append(r"\tablerefs{" + "; ".join(ref_strs) + "}")
    # Footnotes for data-measured vlsr / literature ref
    for tag in sorted(footnotes.keys()):
        lines.append(rf"\tablenotetext{{{tag}}}{{{footnotes[tag]}}}")
    lines.append(r"\tablecomments{Sources are those with $L_\mathrm{bol} \geq 10^4\,L_\odot$ "
                 r"and $d \leq 2$\,kpc drawn from the RMS \citep{Lumsden2013} and HiGAL "
                 r"\citep{Elia2017} catalogs. ALMA project codes were retrieved by querying "
                 r"the ALMA archive for observations within 15\arcsec\ of each position; "
                 r"only those projects whose angular resolution corresponds to a spatial "
                 r"scale finer than 500\,AU at the adopted distance are listed. The "
                 r"$v_\mathrm{LSR}$ column gives the adopted systemic velocity; the "
                 r"$v_\mathrm{LSR}$ reference is `data' when measured from corroborating "
                 r"in-band tracer lines (see per-source footnotes), `RMS' for the Lumsden "
                 r"et al. (2013) tabulated value, or a literature citation otherwise.}")
    lines.append(r"\end{deluxetable}")

    # Companion proposals table — one line per proposal-code per target
    lines.append(r"")
    lines.append(r"\startlongtable")
    lines.append(r"\begin{deluxetable*}{lp{4.5in}}")
    lines.append(r"\tabletypesize{\scriptsize}")
    lines.append(r"\tablecaption{ALMA proposal codes (resolution finer than "
                 r"500\,AU at the source distance) per target.\label{tab:targets_alma}}")
    lines.append(r"\tablehead{\colhead{Source} & "
                 r"\colhead{ALMA proposal codes}}")
    lines.append(r"\startdata")
    rows_with_codes = [r for r in rows if r["codes"]]
    rows_no_codes = [r for r in rows if not r["codes"]]
    for r in rows_with_codes:
        def _fmt_res(arcsec):
            # Use 3 sig figs; "0.05" for 0.045, "0.5" for 0.49, etc.
            if arcsec >= 0.1:
                return f"{arcsec:.2f}"
            return f"{arcsec:.3f}"
        codes_str = ", ".join(
            f"{c}\\,({_fmt_res(res)}\\arcsec)" for c, res in r["per_code"]
        )
        lines.append(f"{r['common']} & {codes_str} \\\\")
    lines.append(r"\enddata")
    lines.append(r"\end{deluxetable*}")

    # Separate table listing sources with no <500-AU ALMA data
    if rows_no_codes:
        lines.append(r"")
        lines.append(r"\startlongtable")
        lines.append(r"\begin{deluxetable}{lllccl}")
        lines.append(r"\tabletypesize{\scriptsize}")
        lines.append(r"\tablecaption{Sample sources with no ALMA observations "
                     r"reaching $<500$\,AU resolution at the source distance. "
                     r"These targets are excluded from the salt-search analysis "
                     r"but listed here for completeness.\label{tab:targets_noalma}}")
        lines.append(r"\tablehead{\colhead{Source} & \colhead{Galactic} & "
                     r"\colhead{R.A.\,Dec.} & "
                     r"\colhead{$d$} & \colhead{$L_\mathrm{bol}$} & "
                     r"\colhead{Reason} \\")
        lines.append(r"& (l\,b) & (J2000) & (kpc) & "
                     r"($10^4\,L_\odot$) & }")
        lines.append(r"\startdata")
        for r in rows_no_codes:
            if r["n_obs_total"] == 0 or r["best_res_au"] is None:
                reason = "No obs"
            else:
                reason = f"best res {r['best_res_au']:.0f}\\,AU"
            lines.append(
                f"{r['common']} & {r['gname']} & \\texttt{{{r['coord']}}} & "
                f"{r['d_str']} & {r['lbol']/1e4:.2f} & {reason} \\\\"
            )
        lines.append(r"\enddata")
        lines.append(r"\end{deluxetable}")

    OUT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {OUT}")
    print(f"  {sum(bool(r['codes']) for r in rows)}/{len(rows)} sources have <500-AU ALMA obs")


if __name__ == "__main__":
    main()
