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
    (231.32083, r"$^{13}$CS 5-4",                "ism"),
    (240.27293, r"CH$_3$OH 5$_{1,4}$-4$_{1,3}$", "ism"),
    (241.79143, r"C$^{34}$S 5-4",                "ism"),
    (244.93556, r"CH$_3$OH 8$_{-1}$-7$_{-1}$",   "ism"),
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

LINE_LISTS = {
    "salt":  NACL_KCL,
    "h2o":   WATER,
    "ism":   ISM,
}

COLOR_MAP = {
    "salt": "red",
    "h2o":  "cyan",
    "ism":  "black",
    "rrl":  "magenta",
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
                 species: tuple = ("salt", "h2o", "ism", "rrl"),
                 arrows: tuple = ("ch3ocho", "ch3oh"),
                 fontsize: int = 9):
    """Overlay rotated text + vertical lines for every line in the active
    species set whose Doppler-shifted observed frequency falls in
    [fmin, fmax].  Vertical arrows mark CH3OCHO + CH3OH XCLASS transitions.

    Returns nothing; caller is expected to have already set xlim/ylim.
    """
    shift = 1.0 - vsys / C_KMS
    used_labels = set()
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
