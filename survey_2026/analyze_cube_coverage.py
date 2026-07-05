"""
Scan downloaded ALMA pipeline FITS cubes (product/*_sci.spw*.cube.I.pbcor.fits)
for each proposal/source and report which salt-relevant lines are covered.

Salt lines to check (rest frequencies):
- NaCl v=0 J=18-17     234.252 GHz
- NaCl v=1 J=18-17     232.510 GHz
- KCl v=0 J=32-31      247.235 GHz
- KCl v=1 J=32-31      245.443 GHz
- H2O 5(1,5)-4(2,2)    325.153 GHz
- H2O 4(1,4)-3(2,1)    232.687 GHz
- H2O 3(1,3)-2(2,0)    183.310 GHz
- H2O 6(1,6)-5(2,3)    22.235  GHz  (maser)
- SiO J=5-4            217.105 GHz
- SiO J=2-1             86.847 GHz

Output table: data/cube_coverage_2026.csv with columns
    Name, proposal_id, MOUS, n_cubes, freq_min_GHz, freq_max_GHz,
    has_NaCl_v0, has_NaCl_v1, has_KCl_v0, has_KCl_v1, has_H2O_232,
    has_H2O_325, has_H2O_183, has_SiO_217, has_SiO_86,
    needs_reimaging, notes

A "needs_reimaging" flag is set if the cube spws don't cover any of our
salt lines (or no cubes were downloaded), meaning the user must
reimage from raw uv data.
"""

import os
import glob
import re
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.io import fits

UVDATA_ROOT = Path("/orange/adamginsburg/salt/survey_2026/uvdata")
SOURCE_TABLE = "data/myso_alma_best_obs.fits"
OUT_CSV = "data/cube_coverage_2026.csv"

# Salt + reference lines (rest GHz, redshift-zero search) — all bands.
# A target gets a "yes" coverage if ANY of these falls inside any spw of any
# downloaded cube, at any band.  Vibrationally-excited NaCl/KCl included.
SALT_LINES = [
    # NaCl ground-state ladder
    ("NaCl_J7",   91.16972),  ("NaCl_J8",  104.18962),
    ("NaCl_J10", 130.22336),  ("NaCl_J14", 182.26628),
    ("NaCl_J17", 221.27269),  ("NaCl_J18", 234.25191),
    ("NaCl_J19", 247.22783),  ("NaCl_J20", 260.20003),
    ("NaCl_J22", 286.13036),  ("NaCl_J25", 324.98864),
    ("NaCl_J27", 350.86819),
    # NaCl vibrationally excited (B6 only, common in hot cores)
    ("NaCl_v1_18-17", 232.50995), ("NaCl_v1_17-16", 219.61494),
    ("NaCl_v2_17-16", 217.97997),
    # KCl ground-state ladder
    ("KCl_J12",   92.10080),  ("KCl_J13",   99.77381),
    ("KCl_J20",  153.42944),  ("KCl_J28",  215.00830),
    ("KCl_J30",  230.32056),  ("KCl_J32",  245.62547),
    ("KCl_J34",  260.91620),  ("KCl_J38",  291.46028),
    ("KCl_J45",  344.78876),
    # KCl vibrational
    ("KClv3_29",  218.57971), ("KClv3_31",  233.60570),
    # H2O across bands
    ("H2O_232",  232.6867),  ("H2O_325",  325.1530),
    ("H2O_183",  183.3101),  ("H2O_336",  336.2280),
    ("H2O_321",  321.2256),  ("H2O_437",  437.3467),
    # Reference SiO
    ("SiO_2-1",   86.84698), ("SiO_5-4",  217.10498),
    ("SiO_8-7",  347.33058),
]
SOURCE_VLSR_KMS = 0.0   # could pull from RMS catalog per source
TOL_GHZ = 0.01          # 10 MHz cushion at the band edges


def cube_freq_range(fits_path):
    """Return (f_min, f_max) in GHz for a 4D ALMA cube FITS header."""
    try:
        hdr = fits.getheader(fits_path)
    except Exception:
        return None, None
    # Spectral axis: usually CTYPE3='FREQ' on ALMA pipeline cubes
    naxis = int(hdr.get("NAXIS", 0))
    for ax in range(1, naxis + 1):
        if hdr.get(f"CTYPE{ax}", "").startswith("FREQ"):
            n = int(hdr.get(f"NAXIS{ax}", 1))
            crval = float(hdr.get(f"CRVAL{ax}", np.nan))   # Hz
            crpix = float(hdr.get(f"CRPIX{ax}", 1))
            cdelt = float(hdr.get(f"CDELT{ax}", np.nan))
            f0 = (crval + (1 - crpix) * cdelt) / 1e9
            f1 = (crval + (n - crpix) * cdelt) / 1e9
            return min(f0, f1), max(f0, f1)
    return None, None


def main():
    obs = Table.read(SOURCE_TABLE)
    rows = []

    for r in obs:
        name = r["Name"].decode() if isinstance(r["Name"], bytes) else r["Name"]
        pid = r["proposal_id"].decode() if isinstance(r["proposal_id"], bytes) else r["proposal_id"]
        src_dir = UVDATA_ROOT / pid / name
        cubes = []
        if src_dir.is_dir():
            cubes = sorted(glob.glob(
                str(src_dir / "**" / "*_sci.spw*.cube.I.pbcor.fits"),
                recursive=True))
            # also scan flat layout (when files are stored alongside .DONE)
            cubes += sorted(glob.glob(
                str(src_dir / "*_sci.spw*.cube.I.pbcor.fits")))
        cubes = sorted(set(cubes))

        if not cubes:
            row = {
                "Name": name, "proposal_id": pid,
                "n_cubes": 0,
                "freq_min_GHz": np.nan, "freq_max_GHz": np.nan,
                "needs_reimaging": True,
                "notes": "no cubes on disk",
            }
            for ln, _ in SALT_LINES:
                row[f"has_{ln}"] = False
            rows.append(row)
            continue

        # Freq ranges
        ranges = []
        for c in cubes:
            f0, f1 = cube_freq_range(c)
            if f0 is not None:
                ranges.append((f0, f1))

        if not ranges:
            row = {
                "Name": name, "proposal_id": pid,
                "n_cubes": len(cubes),
                "freq_min_GHz": np.nan, "freq_max_GHz": np.nan,
                "needs_reimaging": True,
                "notes": "could not parse cube frequencies",
            }
            for ln, _ in SALT_LINES:
                row[f"has_{ln}"] = False
            rows.append(row)
            continue

        f_lo = min(r[0] for r in ranges)
        f_hi = max(r[1] for r in ranges)

        row = {
            "Name": name, "proposal_id": pid,
            "n_cubes": len(cubes),
            "freq_min_GHz": f_lo, "freq_max_GHz": f_hi,
            "needs_reimaging": False,
            "notes": f"{len(ranges)} cubes; {len(set((round(r[0],2), round(r[1],2)) for r in ranges))} unique spws",
        }
        any_line = False
        for ln, nu0 in SALT_LINES:
            covered = any(r[0] - TOL_GHZ <= nu0 <= r[1] + TOL_GHZ for r in ranges)
            row[f"has_{ln}"] = covered
            any_line = any_line or covered
        if not any_line:
            row["needs_reimaging"] = True
            row["notes"] += "; no salt lines in covered spws"
        rows.append(row)

    out = Table(rows=rows)
    out.write(OUT_CSV, format="csv", overwrite=True)
    write_markdown_report(rows)
    print(f"Wrote {OUT_CSV}: {len(out)} sources")
    have = sum(1 for r in rows if not r["needs_reimaging"])
    print(f"  imaged-cube coverage of >=1 salt line : {have}")
    print(f"  needs reimaging                       : {len(rows)-have}")
    for ln, _ in SALT_LINES:
        n = sum(1 for r in rows if r.get(f"has_{ln}"))
        print(f"  has_{ln:9s} : {n}/{len(rows)}")


MD_PATH = "data/cube_coverage_2026.md"


def write_markdown_report(rows):
    rows_sorted = sorted(
        rows,
        key=lambda r: (r["needs_reimaging"], r["proposal_id"], r["Name"]),
    )
    lines = [
        "# Salt-survey ALMA imaged-cube coverage",
        "",
        "Per-source summary of which pipeline-imaged FITS cubes have been",
        "downloaded under `uvdata/<proposal>/<source>/` and which salt-relevant",
        "lines are covered by their spw frequency footprints.",
        "",
        "Salt lines checked (rest GHz):",
        "",
    ] + [f"- {ln}: {nu:.4f}" for ln, nu in SALT_LINES] + [
        "",
        "## Per-source coverage",
        "",
        "| Source | Proposal | n_cubes | freq range (GHz) | covered lines | needs reimaging | notes |",
        "|---|---|---:|---|---|:-:|---|",
    ]
    for r in rows_sorted:
        covered = [ln for ln, _ in SALT_LINES if r.get(f"has_{ln}")]
        f0 = r["freq_min_GHz"]
        f1 = r["freq_max_GHz"]
        frange = (f"{f0:.2f}-{f1:.2f}"
                  if isinstance(f0, float) and not (np.isnan(f0))
                  else "-")
        lines.append(
            f"| {r['Name']} | {r['proposal_id']} | {r['n_cubes']} | "
            f"{frange} | {', '.join(covered) if covered else '-'} | "
            f"{'yes' if r['needs_reimaging'] else 'no'} | {r['notes']} |"
        )

    n_yes = sum(1 for r in rows if r["needs_reimaging"])
    n_no  = len(rows) - n_yes
    lines += [
        "",
        f"**Summary:** {n_no}/{len(rows)} sources have at least one salt line "
        f"covered by an existing pipeline cube; the remaining "
        f"{n_yes}/{len(rows)} need reimaging from raw uv data "
        f"(use `download_alma_data.py --include-asdm` then "
        f"`imaging/run_imaging.py`).",
        "",
    ]
    with open(MD_PATH, "w") as fh:
        fh.write("\n".join(lines))
    print(f"Wrote {MD_PATH}")


if __name__ == "__main__":
    main()
