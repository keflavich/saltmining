"""Shared line-ID style for ALL spectrum plots (lineid_full, spectrum_panels,
kinematic_stack).

Provides:
  LINE_LISTS  -- per-species rest-frame frequency tables
  COLOR_MAP   -- species class -> color
  apply_labels(ax, fmin, fmax, ymax, vsys=0.0, ...) -- overlay rotated text
                labels + vertical lines, Doppler-shifted by vsys.
"""
from __future__ import annotations

import csv
from pathlib import Path

# rest_GHz, label, species_class
NACL_KCL = [
    (217.97715, r"NaCl$_{v=2}$ 17-16", "salt"),
    (219.61300, r"NaCl$_{v=1}$ 17-16", "salt"),
    (220.20640, r"KCl$_{v=2}$ 28-27",  "salt"),
    (221.84200, r"KCl$_{v=1}$ 29-28",  "salt"),
    (228.51760, r"KCl$_{v=2}$ 30-29",  "salt"),
    (230.77580, r"NaCl$_{v=2}$ 18-17", "salt"),
    (232.50782, r"NaCl$_{v=1}$ 18-17", "salt"),
    (234.25090, r"NaCl$_{v=0}$ 18-17", "salt"),
    (243.57049, r"NaCl$_{v=2}$ 19-18", "salt"),
]

WATER = [
    (232.68673, r"H$_2$O 5$_{1,5}$-4$_{2,2}$", "h2o"),
    (232.93657, r"H$_2$O$_{v_2=1}$ 3$_{1,3}$-2$_{2,0}$", "h2o"),
    (234.14050, r"H$_2$O$_{v_2=1}$ 4$_{1,4}$-3$_{2,1}$", "h2o"),
]

# Prominent ISM lines, black labels per user instruction
ISM = [
    (216.11258, r"DCO$^+$ 3-2",                  "ism"),
    (217.10498, r"SiO 5-4",                      "ism"),
    (218.22219, r"H$_2$CO 3$_{0,3}$-2$_{0,2}$",  "ism"),
    (218.47564, r"H$_2$CO 3$_{2,2}$-2$_{2,1}$",  "ism"),  # 322
    (218.76007, r"H$_2$CO 3$_{2,1}$-2$_{2,0}$",  "ism"),  # 321
    (219.56035, r"C$^{18}$O 2-1",                "ism"),
    (220.39870, r"$^{13}$CO 2-1",                "ism"),
    (220.74726, r"CH$_3$CN 12$_3$-11$_3$",       "ism"),
    (220.59442, r"CH$_3$OH 8$_{0}$-7$_{1}$",     "ism"),
    (230.53800, r"$^{12}$CO 2-1",                "ism"),
    (231.06099, r"OCS 19-18",                    "ism"),
    (231.22069, r"$^{13}$CS 5-4",                "ism"),  # was mislabeled at 231.32083
    (240.27293, r"CH$_3$OH 5$_{1,4}$-4$_{1,3}$", "ism"),
    (241.79143, r"C$^{34}$S 5-4",                "ism"),
    (244.93556, r"CH$_3$OH 8$_{-1}$-7$_{-1}$",   "ism"),
    # contaminants confirmed against XCLASS catalogs (2026-07-02); these are
    # DETECTION-GATED in apply_labels (labeled only when the plotted spectrum
    # shows |signal| >= det_nsigma at the shifted frequency)
    (217.23854, r"DCN 3-2",                              "ism"),
    (217.29920, r"CH$_3$OH$_{v_t=1}$ 6$_{1,5}$-7$_{2,6}$", "ism"),
    (217.64286, r"CH$_3$OH$_{v_t=1}$ 15$_{6}$-16$_{5}$",  "ism"),
    (217.88654, r"CH$_3$OH 20$_{1,19}$-20$_{0,20}$",      "ism"),
    (218.68420, r"CH$_3$OCHO$_{v18=1}$ (blend)",          "ism"),
    (218.90336, r"OCS 18-17",                             "ism"),
    (219.17376, r"HC$_3$N$_{v_7=1}$ 24-23 $\ell$=1f",     "ism"),
    (219.65677, r"HNCO 10$_{3,8}$-9$_{3,7}$",             "ism"),
    (219.73719, r"HNCO 10$_{2,8}$-9$_{2,7}$",             "ism"),
    # 231-234 GHz contaminant set (Splatalogue-verified 2026-07-02):
    # NH2CHO J=11-10 K-ladder + CH3OH + HCOOH + SO2. The two low-Eu
    # 'sibling' lines at 218.459 (NH2CHO) and 218.939 (HCOOH) are the
    # falsifiability check: a real formamide/formic detection at 233 GHz
    # REQUIRES these brighter lines to be present where covered.
    (218.45943, r"NH$_2$CHO 10$_{1,9}$-9$_{1,8}$",        "ism"),
    (218.93851, r"HCOOH 8$_{1,8}$-7$_{0,7}$",             "ism"),
    (231.28115, r"CH$_3$OH 10$_{2,9}$-9$_{3,6}$",         "ism"),
    # 231.506: HCOOH 10(1,9)-9(1,8) rejected by the sibling test (HCOOH
    # 218.939, Eu=40K, is absent at 0.3 sigma where covered) -> CH3CHO
    (231.50629, r"CH$_3$CHO 12$_{4,9}$-11$_{4,8}$ (HCOOH?)", "ism"),
    (231.87581, r"U231.876",                              "ism"),
    (232.27390, r"NH$_2$CHO 11$_{2,10}$-10$_{2,9}$",      "ism"),
    (232.78359, r"CH$_3$OH 18$_{3,16}$-17$_{4,13}$",      "ism"),
    (232.94584, r"CH$_3$OH 10$_{-3,8}$-11$_{-2,10}$",     "ism"),
    (233.52830, r"NH$_2$CHO 11$_{6}$-10$_{6}$",           "ism"),
    (233.59500, r"NH$_2$CHO 11$_{5}$-10$_{5}$",           "ism"),
    (233.73516, r"NH$_2$CHO 11$_{4,8}$-10$_{4,7}$",       "ism"),
    (233.74606, r"NH$_2$CHO 11$_{4,7}$-10$_{4,6}$",       "ism"),
    (233.79580, r"CH$_3$OH 18$_{3,15}$-17$_{4,14}$",      "ism"),
    (233.89695, r"NH$_2$CHO 11$_{3,9}$-10$_{3,8}$",       "ism"),
    (234.18705, r"SO$_2$ 28$_{3,25}$-28$_{2,26}$",        "ism"),
    (234.31588, r"NH$_2$CHO 11$_{3,8}$-10$_{3,7}$",       "ism"),
]

# RRLs (sample of useful ones in Band 6 / 7)
RRLS = {
    "H26$\\alpha$": 353.6228,
    "H27$\\alpha$": 316.4156,
    "H28$\\alpha$": 283.5471,
    "H29$\\alpha$": 256.3022,
    "H30$\\alpha$": 231.9009,
    "H31$\\alpha$": 210.5018,
    "H32$\\alpha$": 192.0193,
    "H39$\\beta$":  308.1064,
    "H40$\\beta$":  286.5181,
    "H41$\\beta$":  267.6005,
}

# CH3CN main-isotopologue K-ladder (CDMS, v=0): ALWAYS labeled (all J,K
# in band), per user instruction -- not detection-gated.
CH3CN_LADDER = [
    (202.35551, r"CH$_3$CN 11$_{0}$-10$_{0}$", "ch3cn"),
    (202.35161, r"CH$_3$CN 11$_{1}$-10$_{1}$", "ch3cn"),
    (202.33992, r"CH$_3$CN 11$_{2}$-10$_{2}$", "ch3cn"),
    (202.32044, r"CH$_3$CN 11$_{3}$-10$_{3}$", "ch3cn"),
    (202.29318, r"CH$_3$CN 11$_{4}$-10$_{4}$", "ch3cn"),
    (202.25815, r"CH$_3$CN 11$_{5}$-10$_{5}$", "ch3cn"),
    (202.21537, r"CH$_3$CN 11$_{6}$-10$_{6}$", "ch3cn"),
    (202.16485, r"CH$_3$CN 11$_{7}$-10$_{7}$", "ch3cn"),
    (202.10661, r"CH$_3$CN 11$_{8}$-10$_{8}$", "ch3cn"),
    (202.04068, r"CH$_3$CN 11$_{9}$-10$_{9}$", "ch3cn"),
    (201.96708, r"CH$_3$CN 11$_{10}$-10$_{10}$", "ch3cn"),
    (220.74726, r"CH$_3$CN 12$_{0}$-11$_{0}$", "ch3cn"),
    (220.74301, r"CH$_3$CN 12$_{1}$-11$_{1}$", "ch3cn"),
    (220.73026, r"CH$_3$CN 12$_{2}$-11$_{2}$", "ch3cn"),
    (220.70902, r"CH$_3$CN 12$_{3}$-11$_{3}$", "ch3cn"),
    (220.67929, r"CH$_3$CN 12$_{4}$-11$_{4}$", "ch3cn"),
    (220.64108, r"CH$_3$CN 12$_{5}$-11$_{5}$", "ch3cn"),
    (220.59442, r"CH$_3$CN 12$_{6}$-11$_{6}$", "ch3cn"),
    (220.53932, r"CH$_3$CN 12$_{7}$-11$_{7}$", "ch3cn"),
    (220.47581, r"CH$_3$CN 12$_{8}$-11$_{8}$", "ch3cn"),
    (220.40390, r"CH$_3$CN 12$_{9}$-11$_{9}$", "ch3cn"),
    (220.32363, r"CH$_3$CN 12$_{10}$-11$_{10}$", "ch3cn"),
    (220.23503, r"CH$_3$CN 12$_{11}$-11$_{11}$", "ch3cn"),
    (239.13792, r"CH$_3$CN 13$_{0}$-12$_{0}$", "ch3cn"),
    (239.13331, r"CH$_3$CN 13$_{1}$-12$_{1}$", "ch3cn"),
    (239.11950, r"CH$_3$CN 13$_{2}$-12$_{2}$", "ch3cn"),
    (239.09650, r"CH$_3$CN 13$_{3}$-12$_{3}$", "ch3cn"),
    (239.06430, r"CH$_3$CN 13$_{4}$-12$_{4}$", "ch3cn"),
    (239.02292, r"CH$_3$CN 13$_{5}$-12$_{5}$", "ch3cn"),
    (238.97239, r"CH$_3$CN 13$_{6}$-12$_{6}$", "ch3cn"),
    (238.91272, r"CH$_3$CN 13$_{7}$-12$_{7}$", "ch3cn"),
    (238.84393, r"CH$_3$CN 13$_{8}$-12$_{8}$", "ch3cn"),
    (238.76605, r"CH$_3$CN 13$_{9}$-12$_{9}$", "ch3cn"),
    (238.67912, r"CH$_3$CN 13$_{10}$-12$_{10}$", "ch3cn"),
    (238.58316, r"CH$_3$CN 13$_{11}$-12$_{11}$", "ch3cn"),
    (238.47822, r"CH$_3$CN 13$_{12}$-12$_{12}$", "ch3cn"),
    (257.52738, r"CH$_3$CN 14$_{0}$-13$_{0}$", "ch3cn"),
    (257.52243, r"CH$_3$CN 14$_{1}$-13$_{1}$", "ch3cn"),
    (257.50756, r"CH$_3$CN 14$_{2}$-13$_{2}$", "ch3cn"),
    (257.48279, r"CH$_3$CN 14$_{3}$-13$_{3}$", "ch3cn"),
    (257.44813, r"CH$_3$CN 14$_{4}$-13$_{4}$", "ch3cn"),
    (257.40358, r"CH$_3$CN 14$_{5}$-13$_{5}$", "ch3cn"),
    (257.34918, r"CH$_3$CN 14$_{6}$-13$_{6}$", "ch3cn"),
    (257.28494, r"CH$_3$CN 14$_{7}$-13$_{7}$", "ch3cn"),
    (257.21088, r"CH$_3$CN 14$_{8}$-13$_{8}$", "ch3cn"),
    (257.12704, r"CH$_3$CN 14$_{9}$-13$_{9}$", "ch3cn"),
    (257.03344, r"CH$_3$CN 14$_{10}$-13$_{10}$", "ch3cn"),
    (256.93014, r"CH$_3$CN 14$_{11}$-13$_{11}$", "ch3cn"),
    (256.81716, r"CH$_3$CN 14$_{12}$-13$_{12}$", "ch3cn"),
    (256.69456, r"CH$_3$CN 14$_{13}$-13$_{13}$", "ch3cn"),
    (275.91557, r"CH$_3$CN 15$_{0}$-14$_{0}$", "ch3cn"),
    (275.91026, r"CH$_3$CN 15$_{1}$-14$_{1}$", "ch3cn"),
    (275.89434, r"CH$_3$CN 15$_{2}$-14$_{2}$", "ch3cn"),
    (275.86781, r"CH$_3$CN 15$_{3}$-14$_{3}$", "ch3cn"),
    (275.83068, r"CH$_3$CN 15$_{4}$-14$_{4}$", "ch3cn"),
    (275.78297, r"CH$_3$CN 15$_{5}$-14$_{5}$", "ch3cn"),
    (275.72470, r"CH$_3$CN 15$_{6}$-14$_{6}$", "ch3cn"),
    (275.65589, r"CH$_3$CN 15$_{7}$-14$_{7}$", "ch3cn"),
    (275.57657, r"CH$_3$CN 15$_{8}$-14$_{8}$", "ch3cn"),
    (275.48677, r"CH$_3$CN 15$_{9}$-14$_{9}$", "ch3cn"),
    (275.38653, r"CH$_3$CN 15$_{10}$-14$_{10}$", "ch3cn"),
    (275.27588, r"CH$_3$CN 15$_{11}$-14$_{11}$", "ch3cn"),
    (275.15487, r"CH$_3$CN 15$_{12}$-14$_{12}$", "ch3cn"),
    (275.02356, r"CH$_3$CN 15$_{13}$-14$_{13}$", "ch3cn"),
    (274.88198, r"CH$_3$CN 15$_{14}$-14$_{14}$", "ch3cn"),
    (294.30239, r"CH$_3$CN 16$_{0}$-15$_{0}$", "ch3cn"),
    (294.29673, r"CH$_3$CN 16$_{1}$-15$_{1}$", "ch3cn"),
    (294.27975, r"CH$_3$CN 16$_{2}$-15$_{2}$", "ch3cn"),
    (294.25146, r"CH$_3$CN 16$_{3}$-15$_{3}$", "ch3cn"),
    (294.21187, r"CH$_3$CN 16$_{4}$-15$_{4}$", "ch3cn"),
    (294.16100, r"CH$_3$CN 16$_{5}$-15$_{5}$", "ch3cn"),
    (294.09887, r"CH$_3$CN 16$_{6}$-15$_{6}$", "ch3cn"),
    (294.02550, r"CH$_3$CN 16$_{7}$-15$_{7}$", "ch3cn"),
    (293.94092, r"CH$_3$CN 16$_{8}$-15$_{8}$", "ch3cn"),
    (293.84516, r"CH$_3$CN 16$_{9}$-15$_{9}$", "ch3cn"),
    (293.73828, r"CH$_3$CN 16$_{10}$-15$_{10}$", "ch3cn"),
    (293.62029, r"CH$_3$CN 16$_{11}$-15$_{11}$", "ch3cn"),
    (293.49127, r"CH$_3$CN 16$_{12}$-15$_{12}$", "ch3cn"),
    (293.35125, r"CH$_3$CN 16$_{13}$-15$_{13}$", "ch3cn"),
    (293.20028, r"CH$_3$CN 16$_{14}$-15$_{14}$", "ch3cn"),
    (293.03844, r"CH$_3$CN 16$_{15}$-15$_{15}$", "ch3cn"),
    (312.68774, r"CH$_3$CN 17$_{0}$-16$_{0}$", "ch3cn"),
    (312.68173, r"CH$_3$CN 17$_{1}$-16$_{1}$", "ch3cn"),
    (312.66370, r"CH$_3$CN 17$_{2}$-16$_{2}$", "ch3cn"),
    (312.63365, r"CH$_3$CN 17$_{3}$-16$_{3}$", "ch3cn"),
    (312.59161, r"CH$_3$CN 17$_{4}$-16$_{4}$", "ch3cn"),
    (312.53758, r"CH$_3$CN 17$_{5}$-16$_{5}$", "ch3cn"),
    (312.47158, r"CH$_3$CN 17$_{6}$-16$_{6}$", "ch3cn"),
    (312.39366, r"CH$_3$CN 17$_{7}$-16$_{7}$", "ch3cn"),
    (312.30383, r"CH$_3$CN 17$_{8}$-16$_{8}$", "ch3cn"),
    (312.20213, r"CH$_3$CN 17$_{9}$-16$_{9}$", "ch3cn"),
    (312.08860, r"CH$_3$CN 17$_{10}$-16$_{10}$", "ch3cn"),
    (311.96329, r"CH$_3$CN 17$_{11}$-16$_{11}$", "ch3cn"),
    (311.82625, r"CH$_3$CN 17$_{12}$-16$_{12}$", "ch3cn"),
    (311.67754, r"CH$_3$CN 17$_{13}$-16$_{13}$", "ch3cn"),
    (311.51720, r"CH$_3$CN 17$_{14}$-16$_{14}$", "ch3cn"),
    (311.34530, r"CH$_3$CN 17$_{15}$-16$_{15}$", "ch3cn"),
    (311.16191, r"CH$_3$CN 17$_{16}$-16$_{16}$", "ch3cn"),
    (331.07154, r"CH$_3$CN 18$_{0}$-17$_{0}$", "ch3cn"),
    (331.06518, r"CH$_3$CN 18$_{1}$-17$_{1}$", "ch3cn"),
    (331.04610, r"CH$_3$CN 18$_{2}$-17$_{2}$", "ch3cn"),
    (331.01430, r"CH$_3$CN 18$_{3}$-17$_{3}$", "ch3cn"),
    (330.96979, r"CH$_3$CN 18$_{4}$-17$_{4}$", "ch3cn"),
    (330.91261, r"CH$_3$CN 18$_{5}$-17$_{5}$", "ch3cn"),
    (330.84276, r"CH$_3$CN 18$_{6}$-17$_{6}$", "ch3cn"),
    (330.76028, r"CH$_3$CN 18$_{7}$-17$_{7}$", "ch3cn"),
    (330.66521, r"CH$_3$CN 18$_{8}$-17$_{8}$", "ch3cn"),
    (330.55757, r"CH$_3$CN 18$_{9}$-17$_{9}$", "ch3cn"),
    (330.43741, r"CH$_3$CN 18$_{10}$-17$_{10}$", "ch3cn"),
    (330.30479, r"CH$_3$CN 18$_{11}$-17$_{11}$", "ch3cn"),
    (330.15975, r"CH$_3$CN 18$_{12}$-17$_{12}$", "ch3cn"),
    (330.00234, r"CH$_3$CN 18$_{13}$-17$_{13}$", "ch3cn"),
    (329.83264, r"CH$_3$CN 18$_{14}$-17$_{14}$", "ch3cn"),
    (329.65071, r"CH$_3$CN 18$_{15}$-17$_{15}$", "ch3cn"),
    (329.45661, r"CH$_3$CN 18$_{16}$-17$_{16}$", "ch3cn"),
    (329.25042, r"CH$_3$CN 18$_{17}$-17$_{17}$", "ch3cn"),
    (349.45370, r"CH$_3$CN 19$_{0}$-18$_{0}$", "ch3cn"),
    (349.44699, r"CH$_3$CN 19$_{1}$-18$_{1}$", "ch3cn"),
    (349.42685, r"CH$_3$CN 19$_{2}$-18$_{2}$", "ch3cn"),
    (349.39330, r"CH$_3$CN 19$_{3}$-18$_{3}$", "ch3cn"),
    (349.34634, r"CH$_3$CN 19$_{4}$-18$_{4}$", "ch3cn"),
    (349.28601, r"CH$_3$CN 19$_{5}$-18$_{5}$", "ch3cn"),
    (349.21231, r"CH$_3$CN 19$_{6}$-18$_{6}$", "ch3cn"),
    (349.12529, r"CH$_3$CN 19$_{7}$-18$_{7}$", "ch3cn"),
    (349.02497, r"CH$_3$CN 19$_{8}$-18$_{8}$", "ch3cn"),
    (348.91140, r"CH$_3$CN 19$_{9}$-18$_{9}$", "ch3cn"),
    (348.78462, r"CH$_3$CN 19$_{10}$-18$_{10}$", "ch3cn"),
    (348.64469, r"CH$_3$CN 19$_{11}$-18$_{11}$", "ch3cn"),
    (348.49166, r"CH$_3$CN 19$_{12}$-18$_{12}$", "ch3cn"),
    (348.32558, r"CH$_3$CN 19$_{13}$-18$_{13}$", "ch3cn"),
    (348.14653, r"CH$_3$CN 19$_{14}$-18$_{14}$", "ch3cn"),
    (347.95456, r"CH$_3$CN 19$_{15}$-18$_{15}$", "ch3cn"),
    (347.74977, r"CH$_3$CN 19$_{16}$-18$_{16}$", "ch3cn"),
    (347.53222, r"CH$_3$CN 19$_{17}$-18$_{17}$", "ch3cn"),
    (347.30199, r"CH$_3$CN 19$_{18}$-18$_{18}$", "ch3cn"),
]

# Auto-generated bright XCLASS transitions (build_xclass_labels.py);
# detection-gated like the rest of the ISM class.
try:
    from analysis.xclass_extra_labels import EXTRA_ISM as _EXTRA_ISM
except ImportError:
    try:
        from xclass_extra_labels import EXTRA_ISM as _EXTRA_ISM
    except ImportError:
        _EXTRA_ISM = []
ISM = ISM + list(_EXTRA_ISM)

LINE_LISTS = {
    "salt":  NACL_KCL,
    "h2o":   WATER,
    "ism":   ISM,
    "ch3cn": CH3CN_LADDER,
}

COLOR_MAP = {
    "salt": "red",
    "h2o":  "cyan",
    "ism":  "black",
    "rrl":  "magenta",
    "ch3cn": "darkgreen",
    "ch3ocho": "red",
    "ch3oh":   "blue",
}

C_KMS = 2.99792458e5


def _load_csv_freqs(name):
    """Load CH3OCHO / CH3OH XCLASS transition catalog (data/<name>_transitions.csv)."""
    p = Path(__file__).resolve().parent.parent / f"data/{name}_transitions.csv"
    if not p.exists():
        return []
    out = []
    with p.open() as fh:
        for row in csv.DictReader(fh):
            try:
                f = float(row["rest_GHz"])
            except (KeyError, ValueError, TypeError):
                continue
            out.append(f)
    return out


CH3OCHO_FREQS = _load_csv_freqs("CH3OCHO")
CH3OH_FREQS   = _load_csv_freqs("CH3OH")


def apply_labels(ax, fmin, fmax, ymax,
                 vsys: float = 0.0,
                 species: tuple = ("salt", "h2o", "ism", "rrl", "ch3cn"),
                 arrows: tuple = ("ch3ocho", "ch3oh"),
                 fontsize: int = 9,
                 spec_freq=None, spec_flux=None, det_nsigma: float = 3.0):
    """Overlay rotated text + vertical lines for every line in the active
    species set whose Doppler-shifted observed frequency falls in
    [fmin, fmax].  Vertical arrows mark CH3OCHO + CH3OH XCLASS transitions.

    When `spec_freq`/`spec_flux` (GHz, brightness) are provided, ISM-class
    contaminant labels are DETECTION-GATED: a label is drawn only if the
    spectrum shows |signal| >= det_nsigma * MAD within +-10 km/s of the
    shifted line frequency (emission OR absorption). Search-target classes
    (salt, h2o, rrl) are always labeled.

    Returns nothing; caller is expected to have already set xlim/ylim.
    """
    import numpy as _np
    shift = 1.0 - vsys / C_KMS
    detect = None
    if spec_freq is not None and spec_flux is not None:
        sf = _np.asarray(spec_freq, dtype=float)
        sl = _np.asarray(spec_flux, dtype=float)
        finite = _np.isfinite(sl)
        if finite.sum() > 10:
            med = float(_np.nanmedian(sl[finite]))
            mad = float(1.4826 * _np.nanmedian(_np.abs(sl[finite] - med)))

            def detect(fobs, _sf=sf, _sl=sl, _med=med, _mad=mad):
                if _mad <= 0:
                    return True
                df = fobs * 10.0 / C_KMS
                m = _np.abs(_sf - fobs) <= df
                if not m.any():
                    return False
                return bool(_np.nanmax(_np.abs(_sl[m] - _med)) >=
                            det_nsigma * _mad)
    used_labels = set()
    _ism_placed = []
    for cls in species:
        col = COLOR_MAP.get(cls, "gray")
        if cls == "rrl":
            for lab, rest in RRLS.items():
                fobs = rest * shift
                if not (fmin <= fobs <= fmax):
                    continue
                ax.axvline(fobs, color=col, linestyle=":", linewidth=0.6,
                           alpha=0.7)
                ax.text(fobs, ymax, lab, rotation=90, color=col,
                        fontsize=fontsize, ha="right", va="top", alpha=0.9)
            continue
        for rest, lab, _ in LINE_LISTS.get(cls, []):
            fobs = rest * shift
            if not (fmin <= fobs <= fmax):
                continue
            if cls == "ism" and detect is not None and not detect(fobs):
                continue
            if cls == "ism" and any(abs(fobs - fprev) < 0.008
                                     for fprev in _ism_placed):
                continue  # min 8 MHz label separation to avoid pileup
            if cls == "ism":
                _ism_placed.append(fobs)
            ax.axvline(fobs, color=col, linestyle=":", linewidth=0.6,
                       alpha=0.7)
            ax.text(fobs, ymax, lab, rotation=90, color=col,
                    fontsize=fontsize, ha="right", va="top", alpha=0.9)
            used_labels.add(lab)
    if not arrows:
        return
    ylim_lo, ylim_hi = ax.get_ylim()
    arrow_y = ylim_hi - 0.08 * (ylim_hi - ylim_lo)
    if "ch3ocho" in arrows:
        col = COLOR_MAP["ch3ocho"]
        for rest in CH3OCHO_FREQS:
            fobs = rest * shift
            if fmin <= fobs <= fmax:
                ax.annotate("", xy=(fobs, ylim_hi - 0.18 * (ylim_hi - ylim_lo)),
                            xytext=(fobs, ylim_hi - 0.02 * (ylim_hi - ylim_lo)),
                            arrowprops=dict(arrowstyle="-|>",
                                              color=col, lw=0.8, alpha=0.8),
                            annotation_clip=True)
    if "ch3oh" in arrows:
        col = COLOR_MAP["ch3oh"]
        for rest in CH3OH_FREQS:
            fobs = rest * shift
            if fmin <= fobs <= fmax:
                ax.annotate("", xy=(fobs, ylim_hi - 0.18 * (ylim_hi - ylim_lo)),
                            xytext=(fobs, ylim_hi - 0.02 * (ylim_hi - ylim_lo)),
                            arrowprops=dict(arrowstyle="-|>",
                                              color=col, lw=0.8, alpha=0.8),
                            annotation_clip=True)
