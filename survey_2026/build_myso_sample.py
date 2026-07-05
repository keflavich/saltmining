"""
Build the distance-limited, gas-rich MYSO sample for the salt survey.

Sample definition
-----------------
1. Start from RMS Survey (Lumsden+ 2013) sources with measured kinematic distance
   and RMS-SED-fit bolometric luminosity.
2. Primary type cut: RMS Type in {YSO, HII/YSO, HII region}.
   Pure HII regions are retained ONLY if they show positive "gas-rich" evidence
   (still-embedded hot core candidate behind an HII shell; Orion Source I
   analogue). Evolved stars, PNe, OH/IR stars, Young/old, Other dropped.
3. Distance cut: d < 3.0 kpc (smallest d that yields >=30 ALMA-visible MYSOs).
4. Luminosity cut: L_bol > 1e4 L_sun (conventional MYSO / early B threshold).
5. Gas-rich evidence required for any HII-region-classified source:
      (a) ATLASGAL 870 um counterpart in Urquhart+2014 RMS-ATLASGAL match
          (covers |l|<60 deg, |b|<1.5 deg), OR
      (b) 6.7 GHz Class II methanol maser within 10" in:
          MMB (Caswell+ 2010, VIII/96) +
          Yang+2019 northern |b|<2 deg survey (J/ApJS/241/18) +
          GLOSTAR (Nguyen+ 2022, J/A+A/666/A59).

Known blind spots
-----------------
- ATLASGAL does not cover Cygnus-X, Cassiopeia, Perseus, Gemini, Monoceros,
  Carina, Sco-Cen (|l|>60 deg). HII regions there may be falsely flagged
  gas-poor. Consider JCMT Plane Survey / SCUBA-2 / Bolocam GPS follow-up for
  these candidates.
- The 6.7 GHz maser surveys cover most of the Galactic plane down to |b|<2 deg
  but the non-detection statistics are not homogeneous.

Output
------
data/myso_sample_2026.fits  - full candidate list with provenance columns
data/myso_sample_2026.csv   - human-readable table
"""

import numpy as np
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u


def sexa_coords(t, ra_col="RAJ2000", de_col="DEJ2000"):
    ra = [r.decode() if isinstance(r, bytes) else r for r in t[ra_col]]
    de = [r.decode() if isinstance(r, bytes) else r for r in t[de_col]]
    return SkyCoord(ra, de, unit=(u.hourangle, u.deg))


def byte_col_to_str(col):
    return np.array([x.decode().strip() if isinstance(x, bytes) else x.strip() for x in col])


# -------------------- RMS base catalog --------------------
rms = Table.read("data/rms.fits")
types = byte_col_to_str(rms["Type"])
names = byte_col_to_str(rms["Name"])
d = np.asarray(rms["Dist"], float)
L = np.asarray(rms["Lbol"], float)

DIST_MAX_KPC = 3.0
L_MIN_LSUN = 1e4

cand_mask = (
    np.isin(types, ["YSO", "HII/YSO", "HII region"])
    & np.isfinite(d)
    & np.isfinite(L)
    & (d < DIST_MAX_KPC)
    & (L > L_MIN_LSUN)
)
cand = rms[cand_mask]
cand_c = sexa_coords(cand)
cand_types = types[cand_mask]
cand_names = names[cand_mask]
cand_d = d[cand_mask]
cand_L = L[cand_mask]
print(f"Base candidates (d<{DIST_MAX_KPC} kpc, L>{L_MIN_LSUN:.0e}): {len(cand)}")

# -------------------- ATLASGAL counterpart --------------------
u14 = Table.read("data/urquhart14_rms_atlasgal_match.fits")
u14_rms = set(byte_col_to_str(u14["RMS"]))
has_atlasgal = np.array([n in u14_rms for n in cand_names])
# Pull logLbol & M class (hot core y/n) and offset for records
u14_map = {byte_col_to_str(u14["RMS"])[i]: u14[i] for i in range(len(u14))}
atlas_offset = np.array(
    [float(u14_map[n]["Offset"]) if n in u14_map else np.nan for n in cand_names]
)

# -------------------- 6.7 GHz methanol masers --------------------
mmb = Table.read("data/mmb_caswell2010.fits")
mmb_c = sexa_coords(mmb)

yang = Table.read("data/yang2019_methanol.fits")
yang_c = SkyCoord(
    l=np.asarray(yang["GLON"], float) * u.deg,
    b=np.asarray(yang["GLAT"], float) * u.deg,
    frame="galactic",
).icrs

glo = Table.read("data/glostar_methanol.fits")
glo_c = sexa_coords(glo)

all_maser = SkyCoord(
    np.concatenate([mmb_c.ra.deg, yang_c.ra.deg, glo_c.ra.deg]) * u.deg,
    np.concatenate([mmb_c.dec.deg, yang_c.dec.deg, glo_c.dec.deg]) * u.deg,
)
maser_src = np.concatenate(
    [["MMB"] * len(mmb_c), ["Yang19"] * len(yang_c), ["GLOSTAR"] * len(glo_c)]
)

idx, sep, _ = match_coordinates_sky(cand_c, all_maser)
MASER_MATCH_ARCSEC = 10.0
has_maser = sep < MASER_MATCH_ARCSEC * u.arcsec
maser_catalog = np.where(has_maser, maser_src[idx], "")
maser_sep_arcsec = sep.to(u.arcsec).value

# -------------------- gas-rich / keep logic --------------------
gas_rich = has_atlasgal | has_maser
keep = (cand_types != "HII region") | gas_rich

# ALMA observability
dec_deg = cand_c.dec.deg
alma_vis = dec_deg < 30

# -------------------- assemble output --------------------
gal = cand_c.galactic
out = Table()
out["Name"] = cand_names
out["Type"] = cand_types
out["RAJ2000"] = cand_c.ra.deg
out["DEJ2000"] = cand_c.dec.deg
out["l_deg"] = gal.l.deg
out["b_deg"] = gal.b.deg
out["d_kpc"] = cand_d
out["Lbol_Lsun"] = cand_L
out["has_ATLASGAL"] = has_atlasgal
out["ATLASGAL_offset_arcsec"] = atlas_offset
out["has_6p7GHz_maser"] = has_maser.astype(bool)
out["maser_catalog"] = maser_catalog
out["maser_sep_arcsec"] = maser_sep_arcsec
out["gas_rich"] = gas_rich
out["keep"] = keep
out["ALMA_visible"] = alma_vis

# Sort by distance then Galactic longitude
out.sort(["d_kpc", "l_deg"])

out.write("data/myso_sample_2026.fits", overwrite=True)
out[out["keep"]].write(
    "data/myso_sample_2026_kept.fits", overwrite=True
)
# CSV (drop booleans as 0/1 for readability)
out.write("data/myso_sample_2026.csv", format="csv", overwrite=True)

print(f"\nSample summary (d<{DIST_MAX_KPC} kpc, L>{L_MIN_LSUN:.0e}):")
print(f"  base candidates           : {len(cand)}")
print(f"  with ATLASGAL counterpart : {has_atlasgal.sum()}")
print(f"  with 6.7 GHz maser        : {has_maser.sum()}")
print(f"  gas-rich (ATLASGAL | maser): {gas_rich.sum()}")
print(f"  kept (YSO | HII/YSO | gas-rich HII): {keep.sum()}")
print(f"  kept & ALMA-visible (dec<+30): {(keep & alma_vis).sum()}")
print(f"  northern (dec>+30, not ALMA): {(keep & ~alma_vis).sum()}")
print()
print("Distance breakdown of kept sample:")
for lo, hi in [(0, 1), (1, 1.5), (1.5, 2), (2, 2.5), (2.5, 3)]:
    m = (cand_d >= lo) & (cand_d < hi) & keep
    print(f"  {lo:.1f} - {hi:.1f} kpc : {m.sum()}")

print("\nWrote: data/myso_sample_2026.fits, .csv, _kept.fits")
