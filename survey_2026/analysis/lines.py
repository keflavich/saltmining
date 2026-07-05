"""Salt and reference rest frequencies.

The mainline list is the same set used in `miriam/notebooks/SaltSearchProcess.ipynb`
and `Orion_ALMA_2016.1.00165.S/analysis/lines.py`, restricted here to the lines
we actually stack on or use as velocity templates.

If you need the full Barton & vibrationally-excited tables, import them from
`Orion_ALMA_2016.1.00165.S/analysis/salt_tables.py` (this module re-exports the
loader for convenience).
"""

import sys
from astropy import units as u

# All-band salt lines — NaCl rotational ladder (J=N-(N-1) at ~13.05*N GHz,
# B0 ≈ 6516.7 MHz) and KCl rotational ladder (~7.69*N GHz, B0 ≈ 3845.9 MHz).
# Vibrationally excited transitions are slightly red-shifted; v=0..6 covered.
# Frequencies from Cabezas+ / Cernicharo / Splatalogue / Barton compilations.
# Use load_full_salt_tables() for full v=0..10 ladders.

NACL_LINES = {
    # Band 3 (~84-116 GHz): J=7-6 ground+vib
    "NaCl_v0_7-6":     91.16972 * u.GHz,
    "NaCl_v1_7-6":     90.49994 * u.GHz,
    "NaCl_v2_7-6":     89.83320 * u.GHz,
    "Na37Cl_v0_7-6":   89.22148 * u.GHz,
    # J=8-7
    "NaCl_v0_8-7":    104.18962 * u.GHz,
    "NaCl_v1_8-7":    103.42455 * u.GHz,
    "Na37Cl_v0_8-7":  101.96642 * u.GHz,

    # Band 4 (125-163 GHz): J=10-9, J=11-10
    "NaCl_v0_10-9":   130.22336 * u.GHz,
    "NaCl_v1_10-9":   129.27721 * u.GHz,
    "Na37Cl_v0_10-9": 127.45452 * u.GHz,
    "NaCl_v0_11-10":  143.23752 * u.GHz,
    "NaCl_v1_11-10":  142.19623 * u.GHz,

    # Band 5 (163-211 GHz): J=13-12, J=14-13, J=15-14
    "NaCl_v0_13-12":  169.25966 * u.GHz,
    "NaCl_v0_14-13":  182.26628 * u.GHz,
    "NaCl_v1_14-13":  180.94110 * u.GHz,
    "NaCl_v0_15-14":  195.27077 * u.GHz,
    "NaCl_v1_15-14":  193.85093 * u.GHz,

    # Band 6 (211-275 GHz): J=17-16, J=18-17, J=19-18
    "NaCl_v0_17-16":  221.27269 * u.GHz,
    "NaCl_v1_17-16":  219.61494 * u.GHz,
    "NaCl_v2_17-16":  217.97997 * u.GHz,
    "NaCl_v3_17-16":  216.41000 * u.GHz,
    "NaCl_v4_17-16":  214.74231 * u.GHz,
    "Na37Cl_v0_17-16":216.51500 * u.GHz,
    "Na37Cl_v1_17-16":214.93860 * u.GHz,
    "NaCl_v0_18-17":  234.25191 * u.GHz,
    "NaCl_v1_18-17":  232.50995 * u.GHz,
    "NaCl_v2_18-17":  230.78706 * u.GHz,
    "Na37Cl_v0_18-17":229.24589 * u.GHz,
    "NaCl_v0_19-18":  247.22783 * u.GHz,
    "NaCl_v1_19-18":  245.39495 * u.GHz,
    "NaCl_v0_20-19":  260.20003 * u.GHz,
    "NaCl_v0_21-20":  273.16760 * u.GHz,

    # Band 7 (275-373 GHz): J=22-21..28-27
    "NaCl_v0_22-21":  286.13036 * u.GHz,
    "NaCl_v0_23-22":  299.08818 * u.GHz,
    "NaCl_v0_24-23":  312.04098 * u.GHz,
    "NaCl_v0_25-24":  324.98864 * u.GHz,
    "NaCl_v0_26-25":  337.93108 * u.GHz,
    "NaCl_v0_27-26":  350.86819 * u.GHz,
    "NaCl_v1_27-26":  348.25337 * u.GHz,
    "NaCl_v0_28-27":  363.79988 * u.GHz,

    # Band 8 (385-500 GHz): J=30-29..38-37
    "NaCl_v0_30-29":  389.64639 * u.GHz,
    "NaCl_v0_32-31":  415.46934 * u.GHz,
    "NaCl_v0_34-33":  441.26586 * u.GHz,
    "NaCl_v0_36-35":  467.03396 * u.GHz,
    "NaCl_v0_38-37":  492.77170 * u.GHz,
}

KCL_LINES = {
    # Band 3 (~84-116 GHz): J=11-10..15-14
    "KCl_v0_12-11":    92.10080 * u.GHz,
    "KCl_v0_13-12":    99.77381 * u.GHz,
    "KCl_v0_14-13":   107.44519 * u.GHz,
    "KCl_v0_15-14":   115.11483 * u.GHz,

    # Band 4 (125-163 GHz)
    "KCl_v0_17-16":   130.44897 * u.GHz,
    "KCl_v0_18-17":   138.11267 * u.GHz,
    "KCl_v0_20-19":   153.42944 * u.GHz,

    # Band 5 (163-211 GHz)
    "KCl_v0_22-21":   168.74098 * u.GHz,
    "KCl_v0_24-23":   184.04646 * u.GHz,
    "KCl_v0_26-25":   199.34556 * u.GHz,

    # Band 6 (211-275 GHz)
    "KCl_v0_28-27":   215.00830 * u.GHz,
    "KCl_v0_30-29":   230.32056 * u.GHz,
    "KCl_v0_31-30":   237.97478 * u.GHz,
    "KCl_v0_32-31":   245.62547 * u.GHz,
    "KCl_v0_34-33":   260.91620 * u.GHz,
    "KCl_v0_36-35":   276.19438 * u.GHz,
    # KCl vibrationally excited — these are critical (Orion SrcI shows v=3-6)
    "KClv=3_29-28":   218.57971 * u.GHz,
    "KClv=4_29-28":   217.22891 * u.GHz,
    "KClv=5_29-28":   215.88373 * u.GHz,
    "KClv=6_29-28":   214.54412 * u.GHz,
    "KClv=3_31-30":   233.60570 * u.GHz,
    "KClv=4_31-30":   232.16300 * u.GHz,
    "KClv=5_31-30":   230.72400 * u.GHz,
    "41KCl_29-28":    217.54318 * u.GHz,
    "41KClv=2_31-30": 229.68229 * u.GHz,
    "K37Clv=1_31-30": 229.81883 * u.GHz,

    # Band 7
    "KCl_v0_38-37":   291.46028 * u.GHz,
    "KCl_v0_40-39":   306.71304 * u.GHz,
    "KCl_v0_42-41":   321.95157 * u.GHz,
    "KCl_v0_45-44":   344.78876 * u.GHz,
    "KCl_v0_47-46":   359.99030 * u.GHz,
    "KClv=6_47-46":   346.87489 * u.GHz,
    "KClv=7_47-46":   344.71476 * u.GHz,
}

WATER_LINES = {
    # Band 6
    "H2O_232":       232.6867 * u.GHz,    # v2=1 5(5,0)-6(4,3)
    # Band 7
    "H2O_325":       325.1530 * u.GHz,    # 5(1,5)-4(2,2)
    # Band 5
    "H2O_183":       183.3101 * u.GHz,    # 3(1,3)-2(2,0)
    # Band 7 maser
    "H2O_336":       336.2280 * u.GHz,    # v2=1 5(2,3)-4(3,2)
    "H2O_321":       321.2256 * u.GHz,    # 10(2,9)-9(3,6) maser
    # Band 8
    "H2O_437":       437.3467 * u.GHz,    # 7(5,3)-6(6,0)
    "H2O_439":       439.1508 * u.GHz,    # 6(4,3)-5(5,0)
    "H2O_443":       443.0182 * u.GHz,    # 7(5,3)-7(4,4)
    # Band 1 maser (rarely covered)
    "H2O_22":        22.235  * u.GHz,
}

REFERENCE_LINES = {
    # Band 3
    "SiO_2-1":        86.84698 * u.GHz,
    "C18O_1-0":      109.78217 * u.GHz,
    "13CO_1-0":      110.20137 * u.GHz,
    # Band 4
    "SiO_3-2":       130.26861 * u.GHz,
    # Band 5
    "SiO_4-3":       173.68831 * u.GHz,
    # Band 6
    "SiO_5-4":       217.10498 * u.GHz,
    "SO65-54":       219.94944 * u.GHz,
    "SiSv=0_12-11":  217.81764 * u.GHz,
    "SiSv=0_13-12":  235.96134 * u.GHz,
    "SiSv=0_14-13":  254.10318 * u.GHz,
    "C18O_2-1":      219.56035 * u.GHz,
    "13CO_2-1":      220.39868 * u.GHz,
    "H2CO_303-202":  218.22219 * u.GHz,
    "CH3OH_4-3":     218.44006 * u.GHz,
    "HC3N_24-23":    218.32479 * u.GHz,
    "HNCO_10-9":     219.79827 * u.GHz,
    "CH3CN_12-11":   220.74726 * u.GHz,
    # H2S (ortho 2_2,0-2_1,1)
    "H2S_2_20-2_11": 216.71044 * u.GHz,
    # H2S (para 1_1,0-1_0,1)
    "H2S_1_10-1_01": 168.76276 * u.GHz,
    # PN (J=2-1 thru J=5-4)
    "PN_2-1":         93.97977 * u.GHz,
    "PN_3-2":        140.96776 * u.GHz,
    "PN_4-3":        187.94939 * u.GHz,
    "PN_5-4":        234.93569 * u.GHz,
    # Band 7
    "SiO_8-7":       347.33058 * u.GHz,
    "C18O_3-2":      329.33055 * u.GHz,
    "13CO_3-2":      330.58797 * u.GHz,
    "CO_3-2":        345.79599 * u.GHz,
}

ALL_LINES = {**NACL_LINES, **KCL_LINES, **WATER_LINES, **REFERENCE_LINES}


def lines_in_band(f_lo, f_hi, families=("NaCl", "KCl")):
    """Return ``{name: rest_freq}`` for every line in the given families that
    falls inside ``[f_lo, f_hi]`` (any unit convertible to GHz)."""
    f_lo = u.Quantity(f_lo).to(u.GHz)
    f_hi = u.Quantity(f_hi).to(u.GHz)
    pools = {"NaCl": NACL_LINES, "KCl": KCL_LINES, "H2O": WATER_LINES,
             "ref": REFERENCE_LINES}
    out = {}
    for fam in families:
        for k, v in pools[fam].items():
            if f_lo <= v.to(u.GHz) <= f_hi:
                out[k] = v
    return out


def load_full_salt_tables():
    """Re-export of the Orion-project canonical Barton tables.

    Use these when you want the full v=0..10 ladder for stacking; the curated
    dicts above are only the subset hand-tuned for survey_2026 visibility.
    """
    sys.path.insert(0, "/orange/adamginsburg/salt/Orion_ALMA_2016.1.00165.S/analysis")
    from salt_tables import (salt_tables, salt_table_names, NaCl, NaClv,
                             KCl, K37Cl, Na37Cl, AlCl, AlF, NaF, AlO,
                             AlOH, NaCN, CaS, CaO)
    return {
        "tables": salt_tables, "names": salt_table_names,
        "NaCl": NaCl, "KCl": KCl, "K37Cl": K37Cl, "Na37Cl": Na37Cl,
        "AlCl": AlCl, "AlF": AlF, "NaF": NaF,
        "AlO": AlO, "AlOH": AlOH, "NaCN": NaCN, "CaS": CaS, "CaO": CaO,
    }
