"""
Literature-override supplement for the MYSO sample.

Why this file exists
--------------------
The RMS-based catalog (build_myso_sample.py) is incomplete in two ways:

1. RMS itself omits several famous nearby massive YSO regions:
     - Orion BN/KL / Source I    (l~209, b~-19; outside RMS Galactic-plane focus)
     - NGC 6334 I, I(N)          (l~351, b~+0.7; gap at l=350-352 in RMS)
     - NGC 6357 G353.2           (same gap)
     - Mon R2 IRS3               (RMS associates the wrong source at d=5.9 kpc)
     - S140 IRS1                 (l=107, b=+5.3 not in RMS)
2. RMS has the source but our automated filter rejects it for fixable reasons:
     - NGC 7538 IRS1   - HII region with no ATLASGAL coverage (l>60 deg) so the
                         "gas-rich" flag was False even though it is one of the
                         most-studied massive YSOs/HC HII regions.
     - AFGL 2591       - d=3.30 kpc (just over the 3.0 cut); lit parallax 3.33.
     - DR21 main       - RMS L_bol fit gives 7.8e3 (wrong); lit L>1e5.
     - DR21(OH)        - RMS L_bol fit gives 4.3e2 (wrong); lit L~1e4.
     - Cep A HW2       - RMS L=8.4e3 (just below cut); lit ~1e4.
     - AFGL 490        - L~5e3 (RMS=1.3e3); included as marginal B-star MYSO.

Each row carries the literature reference for d and L_bol so the override is
auditable.  Output is merged with the RMS-derived kept sample to produce
data/myso_sample_2026_full.fits.
"""

import numpy as np
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
import astropy.units as u

# Literature-override entries.
# (Name, alt_name, l_deg, b_deg, d_kpc, Lbol_Lsun, type, dist_ref, L_ref, why_missing)
SUPPLEMENT = [
    ("OrionKL_SrcI",   "Orion BN/KL / Source I", 209.0137, -19.3816, 0.418, 1.0e5, "YSO",
     "Reid+2025 VLBA parallax", "Menten+2007", "RMS lacks; b=-19 outside plane focus"),
    ("NGC6334-I",      "NGC 6334 I (MM1)",       351.4181,  +0.6478, 1.30,  1.6e5, "HII/YSO",
     "Wu+2014 maser parallax", "Sandell 2000", "RMS gap at l=350-352"),
    ("NGC6334-In",     "NGC 6334 I(N)",          351.4096,  +0.6655, 1.30,  1.0e4, "YSO",
     "Wu+2014 maser parallax", "Hunter+2006", "RMS gap at l=350-352"),
    ("NGC6357-G353",   "NGC 6357 G353.2+0.9",    353.1810,  +0.8859, 1.70,  5.0e4, "HII/YSO",
     "Massi+2015",            "Massi+2015", "RMS gap at l=350-353"),
    ("MonR2-IRS3",     "Mon R2 IRS3",            213.7050, -12.6038, 0.83,  6.5e3, "YSO",
     "Schlafly+2014 extinction", "Henning+1992", "RMS associates wrong d=5.9 kpc source"),
    ("S140-IRS1",      "S140 IRS1",              106.7984,  +5.3120, 0.764, 5.0e3, "YSO",
     "Hirota+2008 maser parallax", "Crampton+1974", "Outside RMS; b=+5.3"),
    # in RMS but filter-rejected:
    ("RMS_NGC7538I",   "NGC 7538 IRS1 (=G111.5423+00.7776)", 111.5423,+0.7776, 2.65, 1.5e5, "HII region",
     "Moscadelli+2009 parallax", "RMS SED", "RMS HII; gas_rich=False because l>60 (no ATLASGAL)"),
    ("RMS_AFGL2591",   "AFGL 2591 (=G078.8867+00.7087)",     78.8867, +0.7087, 3.33, 1.9e5, "YSO",
     "Rygl+2012 parallax",      "RMS SED", "d>3.0 kpc cut by 0.33"),
    ("RMS_DR21main",   "DR21 main (=G081.6802+00.5405A)",    81.6802, +0.5405, 1.50, 1.0e5, "HII region",
     "Rygl+2012 parallax",      "Schneider+2010", "RMS L_bol=7.8e3 underestimated"),
    ("RMS_DR21OH",     "DR21(OH) (=G081.7220+00.5699)",      81.7220, +0.5699, 1.50, 1.0e4, "HII region",
     "Rygl+2012 parallax",      "Mangum+1991", "RMS L_bol=4.3e2 underestimated"),
    ("RMS_CepA",       "Cep A HW2 (=G109.8715+02.1156)",    109.8715, +2.1156, 0.70, 1.0e4, "YSO",
     "Moscadelli+2009 parallax","Patel+2005",   "RMS L=8.4e3 just below cut"),
    ("RMS_AFGL490",    "AFGL 490 (=G141.9996+01.8202)",     141.9996, +1.8202, 1.00, 5.0e3, "YSO",
     "Schneider+2011",          "Schreyer+2002","Marginal B-star MYSO; L~5e3"),
]

rows = []
for name, alt, l, b, d, L, typ, dref, Lref, why in SUPPLEMENT:
    c = SkyCoord(l=l*u.deg, b=b*u.deg, frame="galactic").icrs
    rows.append({
        "Name": name,
        "Type": typ,
        "RAJ2000": float(c.ra.deg),
        "DEJ2000": float(c.dec.deg),
        "l_deg": l,
        "b_deg": b,
        "d_kpc": d,
        "Lbol_Lsun": L,
        "has_ATLASGAL": False,
        "ATLASGAL_offset_arcsec": np.nan,
        "has_6p7GHz_maser": False,
        "maser_catalog": "",
        "maser_sep_arcsec": np.nan,
        "gas_rich": True,            # by construction (literature MYSO)
        "keep": True,
        "ALMA_visible": float(c.dec.deg) < 45.0,
        "alt_name": alt,
        "dist_ref": dref,
        "L_ref": Lref,
        "why_added": why,
        "source": "literature_override",
    })
sup = Table(rows)

# Merge with RMS-derived kept sample
base = Table.read("data/myso_sample_2026_kept.fits")
for col in ["alt_name", "dist_ref", "L_ref", "why_added", "source"]:
    if col not in base.colnames:
        if col == "source":
            base[col] = ["RMS"] * len(base)
        else:
            base[col] = [""] * len(base)

# Drop any RMS rows that the supplement supersedes (DR21 main, DR21(OH), AFGL2591,
# NGC7538, AFGL490 are originally NOT in the kept sample, but Cep A and the
# others might exist with stale L; identify by RA/Dec match within 5")
sup_c = SkyCoord(sup["RAJ2000"]*u.deg, sup["DEJ2000"]*u.deg)
base_c = SkyCoord(base["RAJ2000"]*u.deg, base["DEJ2000"]*u.deg)
# remove base rows within 5" of any supplement row
keep_base = np.ones(len(base), bool)
for i, sc in enumerate(sup_c):
    sep = base_c.separation(sc).arcsec
    keep_base &= (sep > 5.0)
n_replaced = (~keep_base).sum()

merged = vstack([base[keep_base], sup])
merged.sort(["d_kpc", "l_deg"])

merged.write("data/myso_sample_2026_full.fits", overwrite=True)
merged.write("data/myso_sample_2026_full.csv", format="csv", overwrite=True)

print(f"RMS kept             : {len(base)}")
print(f"RMS rows replaced    : {n_replaced}")
print(f"Literature additions : {len(sup)}")
print(f"Merged total         : {len(merged)}")
print()
print("Supplement entries:")
for r in sup:
    print(f"  {r['Name']:18s} {r['alt_name']:35s} d={r['d_kpc']:.2f} L={r['Lbol_Lsun']:.1e} - {r['why_added']}")

print("\nWrote: data/myso_sample_2026_full.fits, .csv")
