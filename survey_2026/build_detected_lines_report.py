"""Extended detected-lines report: every (target, proposal, src) with at
least one in-band line at snr >= 3.5, with contaminant candidates listed
within +/- 30 km/s of any suspect 5sigma rest frequency.

Suspect = single-line drives detection OR NaCl v>=2 detected without v=0
corroboration OR H2O at peak_v offset >5 km/s from vsys.

For each suspect, query splatalogue within +/-10 MHz of the rest freq and
report the candidate molecules (Eup, intensity) that could blend.

Output:
  data/detected_lines_report.csv (per-line rows)
  data/detected_lines_report.md  (human-readable per-source markdown)
"""
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u

warnings.filterwarnings("ignore")

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"
OUT_CSV = ROOT / "data/detected_lines_report.csv"
OUT_MD = ROOT / "data/detected_lines_report.md"
OUT_TEX = Path("/orange/adamginsburg/salt/demography_2026/detected_lines_report.tex")

SNR_REPORT = 3.5     # include in the report
SNR_SUSPECT_CHECK = 5.0  # only contaminant-search for this and above

RE_VEXC = re.compile(r"_v(\d+)_")


def vexc(line):
    m = RE_VEXC.search(line)
    return int(m.group(1)) if m else -1


def brightest_id(prop_dir):
    cont = prop_dir / "continuum_sources.csv"
    if not cont.exists() or cont.stat().st_size == 0:
        return None
    try:
        df = pd.read_csv(cont)
    except pd.errors.EmptyDataError:
        return None
    if df.empty or "peak_Jybeam" not in df.columns:
        return None
    return int(df.loc[df["peak_Jybeam"].idxmax(), "id"])


_splat_cache = {}


def splat_query(rest_GHz, tol_MHz=10.0):
    """Return list of (species, rest_GHz, eu_K, log10_aij) within tol_MHz of
    rest_GHz from splatalogue. Cached. Returns [] on error."""
    key = (round(rest_GHz, 4), tol_MHz)
    if key in _splat_cache:
        return _splat_cache[key]
    try:
        from astroquery.splatalogue import Splatalogue
        from astroquery.splatalogue.utils import minimize_table as mt
    except ImportError:
        _splat_cache[key] = []
        return []
    fmin = (rest_GHz - tol_MHz / 1e3) * u.GHz
    fmax = (rest_GHz + tol_MHz / 1e3) * u.GHz
    try:
        t = Splatalogue.query_lines(fmin, fmax,
                                       energy_max=600, energy_type="eu_k",
                                       only_NRAO_recommended=False)
        tb = mt(t)
    except (ConnectionError, OSError, ValueError, KeyError) as e:
        print(f"  splat err: {e}")
        _splat_cache[key] = []
        return []
    out = []
    for r in tb:
        try:
            species = str(r.get("Species") or r.get("species", "?")).strip()
            rest = float(r.get("Freq", r.get("OFreq", np.nan)))
            eu = float(r.get("EU_K", r.get("Eu (K)", np.nan)))
            aij = float(r.get("log10_Aij", r.get("log<sub>10</sub>(Aij)", np.nan)))
        except (TypeError, ValueError):
            continue
        out.append((species, rest, eu, aij))
    _splat_cache[key] = out
    return out


def find_contaminants(line_name, rest_GHz, peak_v_kms, vsys_kms,
                       tol_kms=15.0):
    """Find candidate contaminants near rest_GHz at the observed peak_v.

    Strategy:
      - The peak appears at observed velocity peak_v_kms
      - That observed velocity corresponds to a true rest frequency of
        rest_obs = rest_GHz * (1 - peak_v_kms/c) (radio convention)
      - But for a confusing line at REAL rest frequency rest_other moving at
        vsys, its observed frequency = rest_other * (1 - vsys/c).
        For the observed peak to be at this same observed freq:
            rest_other * (1 - vsys/c) = rest_GHz * (1 - peak_v/c)
        =>  rest_other = rest_GHz * (1 - peak_v/c) / (1 - vsys/c)
      - Search splatalogue near rest_other within tol_kms (converted to MHz)
    """
    c = 299792.458
    rest_other = rest_GHz * (1.0 - peak_v_kms / c) / (1.0 - vsys_kms / c)
    tol_MHz = tol_kms / c * rest_GHz * 1e3
    cands = splat_query(rest_other, tol_MHz=tol_MHz)
    # Exclude self
    cands = [c_ for c_ in cands
             if abs(c_[1] - rest_GHz) > 0.001]
    # Sort by predicted brightness proxy: log10_Aij - Eu/100 (rough)
    def score(t):
        sp, r, eu, aij = t
        if not np.isfinite(aij): aij = -6
        if not np.isfinite(eu): eu = 500
        return -(aij - eu / 100.0)  # smaller = better (Aij large, Eu small)
    return sorted(cands, key=score)[:5]


def vsys_for(target):
    """Same lookup as the kinematic-stack runner."""
    import json
    p = ROOT / "data/vlsr_from_data.json"
    if p.exists():
        try:
            d = json.loads(p.read_text())
            v = d.get(target, {}).get("v_LSR_kms")
            if v is not None: return float(v)
        except (json.JSONDecodeError, KeyError): pass
    p = ROOT / "data/vlsr_from_literature.json"
    if p.exists():
        try:
            d = json.loads(p.read_text())
            v = d.get(target, {}).get("v_LSR_kms")
            if v is not None: return float(v)
        except (json.JSONDecodeError, KeyError): pass
    manual = {"MonR2-IRS3": 10.0, "MonR2-IRS2": 10.0,
              "NGC6334I": -7.0, "NGC6334IN": -3.5,
              "Orion_SrcI": 5.0, "Orion-SrcI": 5.0,
              "G326.6618+00.5207": -39.6, "G345.5043+00.3480": -17.0}
    return manual.get(target, 0.0)


def main():
    rows = []
    for prop_dir in sorted(ANALYSIS.glob("*/2*")):
        meas = prop_dir / "line_measurements.csv"
        if not meas.exists() or meas.stat().st_size == 0:
            continue
        try:
            df = pd.read_csv(meas)
        except pd.errors.EmptyDataError:
            continue
        if df.empty or "snr" not in df.columns:
            continue
        bid = brightest_id(prop_dir)
        if bid is None:
            continue
        sub = df[df["source"] == bid].copy()
        if "rest_GHz" not in sub.columns:
            continue
        hits = sub[sub["snr"] >= SNR_REPORT]
        if hits.empty:
            continue
        target = prop_dir.parent.name
        proposal = prop_dir.name
        vsys = vsys_for(target)
        for _, r in hits.iterrows():
            line = str(r["line"])
            snr = float(r["snr"])
            peak_v = float(r.get("peak_v", np.nan))
            rest = float(r.get("rest_GHz", np.nan))
            group = str(r.get("group", ""))
            v_state = vexc(line)
            # Suspect criteria
            suspects = []
            if snr >= SNR_SUSPECT_CHECK and v_state >= 2 and group in ("NaCl", "KCl"):
                suspects.append("vexc-only (v>=2)")
            if group == "H2O" and np.isfinite(peak_v) and \
                    abs(peak_v - vsys) > 5:
                suspects.append(f"peak_v offset {peak_v-vsys:+.1f} km/s from vsys")
            # Contaminants
            contaminants = []
            if snr >= SNR_SUSPECT_CHECK and np.isfinite(rest) and np.isfinite(peak_v):
                cands = find_contaminants(line, rest, peak_v, vsys)
                for sp, r_ghz, eu, aij in cands:
                    contaminants.append(f"{sp} ({r_ghz:.4f} GHz, Eu={eu:.0f}K)")
            rows.append({
                "target": target, "proposal": proposal, "src_id": bid,
                "line": line, "group": group, "rest_GHz": rest,
                "snr": snr, "peak_v": peak_v, "v_state": v_state,
                "is_5sig": snr >= SNR_SUSPECT_CHECK,
                "suspect_flags": "; ".join(suspects),
                "contaminant_candidates": "; ".join(contaminants[:3]),
            })
    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV} ({len(out)} rows)")

    # Markdown report
    md = ["# Detected lines report",
          f"\nReporting all (target, proposal, src) with snr >= {SNR_REPORT}σ "
          f"at the brightest mm continuum source. Suspect flags + likely "
          f"contaminants listed for >= {SNR_SUSPECT_CHECK}σ detections.\n"]
    for (tgt, prop, sid), g in out.groupby(["target", "proposal", "src_id"]):
        md.append(f"\n## {tgt}  {prop}  src{int(sid):02d}\n")
        g_sorted = g.sort_values("snr", ascending=False)
        md.append("| Line | Group | snr | peak_v | suspect | likely contaminants |")
        md.append("|---|---|---|---|---|---|")
        for _, r in g_sorted.iterrows():
            cont = r["contaminant_candidates"] or ""
            sus = r["suspect_flags"] or ""
            star = " ★" if r["is_5sig"] else ""
            md.append(f"| {r['line']} | {r['group']} | "
                      f"{r['snr']:.1f}σ{star} | {r['peak_v']:+.1f} | "
                      f"{sus} | {cont} |")
    OUT_MD.write_text("\n".join(md))
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main()
