"""Download specific MOUS pipeline tarballs for salt-line targets that
were imaged under different proposals than appear in our myso_alma_best_obs
table.

Targets:
  GGD 27 / HH 80-81 / IRAS 18162-2048 (G010.8411-02.5919)
    - 2016.1.01036.S MOUSes uid://A001/X88e/Xfe and uid://A001/X88e/X100
    - 2019.2.00120.S MOUS uid://A001/X14ed/X1ea
  IRAS 07299-1651 / RCW 7 (G232.6207+00.9959)
    - 2015.1.01454.S MOUSes uid://A001/X2fb/Xe5b and uid://A001/X2fb/Xe5d
    - 2016.1.00125.S MOUS uid://A001/X87d/X92f

Pulls only ``*_sci*cube*pbcor.fits`` URLs to keep size down.
"""

import os
import sys
from pathlib import Path

from astroquery.alma import Alma

UV = Path("/orange/adamginsburg/salt/survey_2026/uvdata")

TARGETS = [
    # GGD 27 / HH 80-81 (high-res any-band)
    ("2016.1.01036.S", "uid://A001/X88e/Xfe",  "G010.8411-02.5919"),  # B6 0.04"
    ("2017.1.01079.S", "uid://A001/X1288/X5fd","G010.8411-02.5919"),  # B3 0.14" NaCl 104
    ("2015.1.00480.S", "uid://A001/X2d8/X436", "G010.8411-02.5919"),  # B7 0.15" NaCl 312
    ("2017.1.00101.S", "uid://A001/X12a3/X1a7","G010.8411-02.5919"),  # B6 0.23" KCl 245
    # IRAS 07299 (already partially downloaded)
    ("2015.1.01454.S", "uid://A001/X2fb/Xe5b", "G232.6207+00.9959"),
    ("2015.1.01454.S", "uid://A001/X2fb/Xe5d", "G232.6207+00.9959"),
    ("2016.1.00125.S", "uid://A001/X87d/X92f", "G232.6207+00.9959"),
]


def fetch(pid, mous, srcname):
    tgt = UV / pid / srcname
    tgt.mkdir(parents=True, exist_ok=True)
    print(f"\n=== {pid} / {mous} -> {tgt} ===")
    alma = Alma()
    try:
        urls = alma.get_data_info(mous, expand_tarfiles=True)
    except Exception as e:
        print(f"  get_data_info failed: {e}; trying expanded=False")
        urls = alma.get_data_info(mous, expand_tarfiles=False)
    # Filter to science-target spw cube products.  Be permissive on naming:
    # accept any "*cube*pbcor.fits" (covers ari_l, manual, *_sci.spw*.cube.*).
    keep = []
    for row in urls:
        url = row["access_url"] if "access_url" in row.colnames else row["URL"]
        if "cube" in url and url.endswith(".fits"):
            keep.append(url)
    print(f"  matched {len(keep)} science-cube URLs (expanded {len(urls)} total)")
    for u in keep[:5]:
        print(f"    {u}")
    if not keep:
        print("  (no science-cube URLs found; falling back to whole product tar)")
        keep = [r["access_url"] if "access_url" in urls.colnames else r["URL"]
                for r in urls if "_001_of_001.tar" in (
                    r["access_url"] if "access_url" in urls.colnames else r["URL"])]
        for u in keep[:3]:
            print(f"    fallback: {u}")
    if not keep:
        print("  nothing to fetch")
        return
    Alma.cache_location = str(tgt)
    try:
        alma.download_files(keep, savedir=str(tgt))
    except Exception as e:
        print(f"  download failed: {e}")


def main():
    for pid, mous, src in TARGETS:
        fetch(pid, mous, src)


if __name__ == "__main__":
    main()
