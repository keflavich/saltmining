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
    "NGC6334I":            r"\citet{Chibueze2014}",
    "NGC6334IN":           r"\citet{Chibueze2014}",
    "MonR2-IRS3":          r"\citet{Herbst1976}",
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
    "OrionSrcI":          r"\citet{Reid2007}",
    "Orion_SrcI":         r"\citet{Reid2007}",
    "NGC6334I":           r"\citet{Sandell2000}",
    "NGC6334IN":          r"\citet{Sandell2000}",
    "MonR2-IRS3":         r"\citet{Henning1992}",
    "G268.4222-00.8490":  r"\citet{Getman2019}",
}


def query_alma_lt500au(coord, dist_kpc, radius=15 * u.arcsec):
    """Return sorted unique project codes with spatial_resolution * d_pc < 500 AU."""
    try:
        q = Alma.query_region(coord, radius=radius)
    except (ConnectionError, ValueError) as e:
        print(f"  query error: {e}")
        return []
    if q is None or len(q) == 0:
        return []
    sr_col = "spatial_resolution" if "spatial_resolution" in q.colnames else "s_resolution"
    pc_col = "proposal_id" if "proposal_id" in q.colnames else "project_code"
    sr = np.asarray(q[sr_col], dtype=float)
    pc = np.asarray(q[pc_col], dtype=str)
    res_au = sr * dist_kpc * 1000.0
    mask = np.isfinite(sr) & (sr > 0) & (res_au < 500.0)
    return sorted(set(pc[mask]))


def format_proposals(codes, per_line=1):
    """Wrap a list of project codes onto multiple lines for table cell."""
    if not codes:
        return r"\nodata"
    chunks = [codes[i:i + per_line] for i in range(0, len(codes), per_line)]
    lines = [", ".join(c) for c in chunks]
    return r" \newline ".join(lines)


def fmt_coord(ra_deg, dec_deg):
    c = ac.SkyCoord(ra_deg * u.deg, dec_deg * u.deg)
    ra = c.ra.to_string(unit=u.hourangle, sep=":", pad=True, precision=1)
    dec = c.dec.to_string(unit=u.deg, sep=":", pad=True, precision=0, alwayssign=True)
    return ra, dec


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
        codes = query_alma_lt500au(coord, r.dist_kpc)
        ra, dec = fmt_coord(r.ra_deg, r.dec_deg)
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
            ra=ra, dec=dec,
            d=r.dist_kpc, lbol=r.lbol_lsun,
            dref=dref_n, lref=lref_n,
            v_lsr=v_str, vref=vref_n,
            codes=codes,
        ))

    pd.DataFrame([{**r, "codes": "; ".join(r["codes"])} for r in rows]).to_csv(
        ROOT / "data/target_table_alma.csv", index=False
    )

    lines = []
    # Compact portrait version of the targets table. ALMA proposals are
    # moved to a separate table; this one lists meta + vLSR only.
    lines.append(r"\startlongtable")
    lines.append(r"\begin{deluxetable}{lccccccccc}")
    lines.append(r"\tabletypesize{\scriptsize}")
    lines.append(r"\tablecaption{Target source list. ALMA project codes are listed when "
                 r"the corresponding observations achieved a spatial resolution finer than "
                 r"500\,AU at the source distance.\label{tab:targets}}")
    lines.append(r"\tablehead{")
    lines.append(r"\colhead{Source} & \colhead{R.A.} & \colhead{Dec.} & "
                 r"\colhead{$d$} & \colhead{$d^{*}$} & "
                 r"\colhead{$L_\mathrm{bol}$} & \colhead{$L^{*}$} & "
                 r"\colhead{$v_\mathrm{LSR}$} & \colhead{$v^{*}$} & "
                 r"\colhead{$N_\mathrm{ALMA}$} \\")
    lines.append(r"& (J2000) & (J2000) & (kpc) & ref. & ($10^4\,L_\odot$) "
                 r"& ref. & (km\,s$^{-1}$) & ref. & ($<500$\,AU)}")
    lines.append(r"\startdata")
    for r in rows:
        n_alma = len(r["codes"])
        lines.append(
            f"{r['name']} & {r['ra']} & {r['dec']} & "
            f"{r['d']:.2f} & {r['dref']} & "
            f"{r['lbol']/1e4:.2f} & {r['lref']} & "
            f"{r['v_lsr']} & {r['vref']} & "
            f"{n_alma if n_alma > 0 else r'\\nodata'} \\\\"
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
    lines.append(r"\begin{deluxetable}{lp{3.5in}}")
    lines.append(r"\tabletypesize{\scriptsize}")
    lines.append(r"\tablecaption{ALMA proposal codes (resolution finer than "
                 r"500\,AU at the source distance) per target.\label{tab:targets_alma}}")
    lines.append(r"\tablehead{\colhead{Source} & "
                 r"\colhead{ALMA proposal codes}}")
    lines.append(r"\startdata")
    for r in rows:
        if not r["codes"]:
            lines.append(f"{r['name']} & \\nodata \\\\")
        else:
            codes_str = ", ".join(r["codes"])
            lines.append(f"{r['name']} & {codes_str} \\\\")
    lines.append(r"\enddata")
    lines.append(r"\end{deluxetable}")

    OUT.write_text("\n".join(lines) + "\n")
    print(f"\nwrote {OUT}")
    print(f"  {sum(bool(r['codes']) for r in rows)}/{len(rows)} sources have <500-AU ALMA obs")


if __name__ == "__main__":
    main()
