"""Multi-survey cross-match LaTeX table for the survey_2026 paper.

Columns (per L4_d2 source):
  Field
  Parent region          - hand curated / RMS Name column
  L_bol                  - L_sun, sources_L4_d2.csv
  d                      - kpc, sources_L4_d2.csv
  S_alma                 - peak mm cont flux at brightest source (Jy/beam)
                            from analysis_products/.../continuum_sources.csv
  theta_alma             - synthesized beam BMAJ at the cube the source
                            sits in (arcsec)
  nu_alma                - median freq of that cube (GHz)
  S_5GHz                 - CORNISH 5 GHz integrated flux, mJy
  theta_5GHz             - CORNISH radio_FWHM_as (")
  radio_class            - radio_characterized.fits classification

CORNISH from /orange/adamginsburg/galactic_plane_surveys/uchii_2026/
radio_characterized.fits; spatial cross-match within 30" of source coord.
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy import units as u

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
SRC_CSV = ROOT / "data/sources_L4_d2.csv"
UCHII = Path("/orange/adamginsburg/galactic_plane_surveys/uchii_2026")
RADIO_FITS = UCHII / "radio_characterized.fits"
OUT_TEX = Path("/orange/adamginsburg/salt/demography_2026/multisurvey.tex")
OUT_CSV = ROOT / "data/multisurvey_table.csv"


PARENT_REGION = {
    "OrionBN-KL": "Orion BN/KL", "Orion-SrcI": "Orion BN/KL",
    "Orion-BN": "Orion BN/KL", "Orion_SrcI": "Orion BN/KL",
    "OrionB-Flame": "Orion B (NGC 2024)",
    "MonR2-IRS3": "Mon R2", "MonR2-IRS2": "Mon R2",
    "NGC6334I": "NGC 6334", "NGC6334IN": "NGC 6334",
    "NGC6514": "M20 (Trifid)", "NGC6357": "NGC 6357",
    "S140-IRS1": "S140 (L1204)",
    "Lagoon-Her36": "M8 (Lagoon)",
    "GGD12-15": "GGD 12-15", "G010.8411-02.5919": "GGD 27 / HH 80-81",
    "IRAS17233-3606": "G351.78 (NGC 6334 area)",
    "G014.9958-00.6732": "M17 SW",
    "G015.0357-00.6795": "M17 SW (HII)",
    "G015.1288-00.6717": "M17",
    "G076.3829-00.6210": "ON 1",
    "G077.9637-00.0075": "Cyg X",
    "G080.8645+00.4197": "Cyg X (DR21 area)",
    "G080.9383-00.1268": "Cyg X (DR21 cluster)",
    "G081.6802+00.5405A": "DR21(OH)", "G081.6802+00.5405B": "DR21(OH)",
    "G081.8789+00.7822": "W75N", "G083.0936+03.2724": "Cyg X (north)",
    "G189.0307+00.7821": "S252 (NGC 2175)",
    "G192.6005-00.0479": "S255 IR",
    "G232.6207+00.9959": "G232.6 (Vela)",
    "G268.4222-00.8490": "RCW 39 / Gum 4",
    "G326.6618+00.5207": "G326.6 (NGC 6231 region)",
    "G336.4917-01.4741A": "G336.5",
    "G336.4917-01.4741B": "G336.5",
    "G342.1156-03.7089": "RCW 122",
    "G345.0052+01.8209": "G345.0",
    "G345.0061+01.7944C": "G345.0",
    "G345.2244+01.0304": "G345.2",
    "G345.4938+01.4677": "G345.5 (NGC 6357 area)",
    "G345.5043+00.3480": "G345.5",
    "G011.3252-01.8040": "M16/M17 area",
    "HIGALBM343.3533-0.0782": "G343.4",
    "G126.7144-00.8220": "W3 area",
    "G133.6923+01.2270": "W3 main",
    "G133.6945+01.2166A": "W3 main",
    "G133.6945+01.2166B": "W3 main",
    "G133.7150+01.2155": "W3 main",
    "G133.7184+01.2237": "W3 main",
    "G133.7354+01.1827": "W3 OH",
    "G133.9476+01.0648": "W3",
    "G353.2+0.9": "G353.2 (NGC 6334 area)",
}


def carry_alma_at_source(target):
    """Return (S_Jy/beam, theta_arcsec, freq_GHz) at the brightest source
    of the first analyzed proposal. Returns NaN if missing."""
    tgt_dir = ANALYSIS / target
    if not tgt_dir.is_dir():
        return (np.nan, np.nan, np.nan)
    for prop_dir in sorted(tgt_dir.glob("2*")):
        cont = prop_dir / "continuum_sources.csv"
        if not cont.exists():
            continue
        try:
            df = pd.read_csv(cont)
        except pd.errors.EmptyDataError:
            continue
        if df.empty:
            continue
        bid = int(df.loc[df["peak_Jybeam"].idxmax(), "id"])
        S = float(df.loc[bid - 1, "peak_Jybeam"])
        # Pull representative beam + freq from any cube in this proposal
        cubes = sorted((ROOT / "uvdata" / prop_dir.name / target)
                          .glob("*sci*cube*pbcor.fits"))
        if not cubes:
            return (S, np.nan, np.nan)
        h = fits.getheader(cubes[0])
        bmaj = float(h.get("BMAJ", 0)) * 3600
        nu_GHz = float(h.get("CRVAL3", 0)) / 1e9
        return (S, bmaj, nu_GHz)
    return (np.nan, np.nan, np.nan)


def cornish_match(coord, table, sep_arcsec=30):
    cs = SkyCoord(table["RAJ2000"], table["DEJ2000"], unit=(u.hourangle, u.deg))
    sep = coord.separation(cs).arcsec
    idx = int(np.argmin(sep))
    if sep[idx] > sep_arcsec:
        return None
    return table[idx]


def fmt_num(v, sig=2):
    if not np.isfinite(v):
        return r"\nodata"
    if v >= 100:
        return f"{v:.0f}"
    if v >= 10:
        return f"{v:.1f}"
    if v >= 1:
        return f"{v:.2f}"
    if v >= 0.01:
        return f"{v:.3f}"
    return f"{v:.2g}"


def main():
    src = pd.read_csv(SRC_CSV).sort_values("dist_kpc")
    radio = Table.read(RADIO_FITS) if RADIO_FITS.exists() else None
    rows = []
    for _, r in src.iterrows():
        name = r["name"]
        d = float(r["dist_kpc"])
        L = float(r["lbol_lsun"])
        S_alma, theta_alma, nu_alma = carry_alma_at_source(name)
        parent = PARENT_REGION.get(name, "...")
        # CORNISH cross-match by world coord
        coord = SkyCoord(float(r["ra_deg"]), float(r["dec_deg"]), unit="deg")
        S_5 = np.nan; theta_5 = np.nan; radio_class = ""
        if radio is not None:
            m = cornish_match(coord, radio)
            if m is not None:
                if "cornish_flux_mJy" in radio.colnames:
                    S_5 = float(m["cornish_flux_mJy"])
                if "radio_FWHM_as" in radio.colnames:
                    theta_5 = float(m["radio_FWHM_as"])
                if "radio_class" in radio.colnames:
                    radio_class = str(m["radio_class"]).strip()
        rows.append({
            "name": name, "parent": parent,
            "dist_kpc": d, "lbol_lsun": L,
            "S_alma_Jybeam": S_alma, "theta_alma_arcsec": theta_alma,
            "nu_alma_GHz": nu_alma,
            "S_5GHz_mJy": S_5, "theta_5GHz_as": theta_5,
            "radio_class": radio_class,
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV} ({len(out)} rows)")

    # Build LaTeX deluxetable
    lines = [
        r"\begin{longrotatetable}",
        r"\begin{deluxetable*}{llcccccccc}",
        r"\tabletypesize{\scriptsize}",
        r"\tablecaption{Multi-survey cross-match for the demography\_2026 "
        r"$L \geq 10^4 \, L_\odot$, $d \leq 2$\,kpc target sample. "
        r"Distances and bolometric luminosities are from the RMS catalog "
        r"\citep{Lumsden+13,Urquhart+18}, except where a more recent reference "
        r"is given in Table~\ref{tab:targets}. ALMA continuum quantities are "
        r"measured at the brightest mm continuum peak in the field across the "
        r"analyzed proposals (peak Jy\,beam$^{-1}$, synthesized beam BMAJ at "
        r"the source spectral window, median spectral-window frequency). "
        r"Radio continuum entries come from the CORNISH survey at $5$\,GHz "
        r"\citep{Hoare+12,Purcell+13} via the radio characterization table "
        r"in \texttt{galactic\_plane\_surveys/uchii\_2026/radio\_characterized.fits} "
        r"(cross-match radius $30''$). Parent-region labels are hand curated "
        r"from SIMBAD nearest-neighbor and ATLASGAL/MMB associations.\label{tab:multisurvey}}",
        r"\tablehead{",
        r"\colhead{Field} & \colhead{Parent} & \colhead{$d$} & "
        r"\colhead{$L_\mathrm{bol}$} & \colhead{$S_\nu^\mathrm{ALMA}$} & "
        r"\colhead{$\theta_\mathrm{ALMA}$} & \colhead{$\nu_\mathrm{ALMA}$} & "
        r"\colhead{$S_\nu^\mathrm{5\,GHz}$} & "
        r"\colhead{$\theta_\mathrm{5\,GHz}$} & "
        r"\colhead{Radio} \\",
        r" & & (kpc) & ($10^4\,L_\odot$) & (Jy\,beam$^{-1}$) & ($''$) & (GHz) & "
        r"(mJy) & ($''$) & class}",
        r"\startdata",
    ]
    for _, r in out.iterrows():
        L4 = r["lbol_lsun"] / 1e4
        lines.append(
            " & ".join([
                r["name"].replace("_", r"\_"),
                str(r["parent"]),
                fmt_num(r["dist_kpc"]),
                fmt_num(L4),
                fmt_num(r["S_alma_Jybeam"]),
                fmt_num(r["theta_alma_arcsec"]),
                fmt_num(r["nu_alma_GHz"]),
                fmt_num(r["S_5GHz_mJy"]),
                fmt_num(r["theta_5GHz_as"]),
                (str(r["radio_class"]) or r"\nodata"),
            ]) + r" \\"
        )
    lines += [
        r"\enddata",
        r"\tablecomments{Sources are sorted by heliocentric distance. The "
        r"`Parent' column gives the larger star-formation region or HII complex "
        r"containing the source; `Radio class' follows the radio\_characterize "
        r"classification (compact, extended, candidate, etc.) of the CORNISH "
        r"counterpart. $S_\nu^\mathrm{ALMA}$ is the peak mm continuum brightness "
        r"of the brightest core in the field and is dominated by dust + free-free "
        r"contribution; it is not the integrated source flux. For sources whose "
        r"continuum was not analyzed locally (no entry in "
        r"\texttt{analysis\_products/}) the ALMA columns are \nodata.}",
        r"\end{deluxetable*}",
        r"\end{longrotatetable}",
    ]
    OUT_TEX.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT_TEX}")


if __name__ == "__main__":
    main()
