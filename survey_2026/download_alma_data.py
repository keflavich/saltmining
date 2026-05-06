"""
Download ALMA archive data for the salt-survey MYSO targets that we don't
yet have on disk.

Inputs
------
- data/myso_alma_best_obs.fits : (Name, RA/Dec, proposal_id, ...)
- data/myso_alma_disk_inventory.csv : produced by data_disk_inventory.py

Strategy
--------
For each (proposal_id, source) pair flagged on_disk=no, query the ALMA
archive for that *position + project* combination, take the row with the
smallest spatial_resolution (matching the "best" observation we previously
identified), expand the data package URLs (including tarfiles), filter to
science MS files (uid___A002_*.ms / *.ms.tar / *_split_cal*.tgz / *_calibrated*.ms),
and download with Alma.download_and_extract_files.

Storage layout (mirrors dihca2's "original_uvdata" pattern):
  /orange/adamginsburg/salt/survey_2026/uvdata/<proposal_id>/<Name>/<files>

A per-target marker file <name>.DONE is written when the download succeeds.

Usage
-----
  # Authenticate once interactively
  python -c "from astroquery.alma import Alma; Alma().login('keflavich')"

  # Dry run (lists what would be fetched, no files downloaded)
  python download_alma_data.py --dry-run

  # Real download (default downloads only science MS / split-cal files)
  python download_alma_data.py

  # Restrict to specific proposals
  python download_alma_data.py --proposals 2022.1.01344.S 2023.1.01346.S

  # Also fetch product .pbcor.fits cubes (huge; off by default)
  python download_alma_data.py --include-products
"""

import os
import re
import sys
import argparse
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

UVDATA_ROOT = Path("/orange/adamginsburg/salt/survey_2026/uvdata")
INVENTORY_CSV = "data/myso_alma_disk_inventory.csv"
BEST_OBS_FITS = "data/myso_alma_best_obs.fits"

# Default: auxiliary tarball (small, contains scriptForPI / restore
# scripts + caltables) PLUS individual science-target imaged products
# (`*_sci.spw*.pbcor.fits` and `*_sci.spw*.mfs.*pbcor.fits`).  We deliberately
# skip the full project tarball (`_001_of_001.tar`), which bundles every
# imaged calibrator + every spw and is hundreds of GB for big projects.
# This default keeps salt-relevant cubes while pruning calibrator imaging.
MS_PATTERNS = [
    re.compile(r"_auxiliary\.tar$"),
    re.compile(r"_sci\.spw[\d_]+\.cube\.I\.pbcor\.fits$"),
    re.compile(r"_sci\.spw[\d_]+\.mfs\.I\.pbcor\.fits$"),
    re.compile(r"_sci\.spw[\d_]+\.mfs\.I\.tt0?\.pbcor\.fits$"),
    re.compile(r"_sci\.spw[\d_]+\.mfs\.IQUV\.manual\.pbcor\.fits$"),
]
# With --include-project-tar: also fetch the full Member OUS product tarball
# (HUGE; bundles every spw + every calibrator imaging).
PROJECT_TAR_PATTERNS = [
    re.compile(r"_001_of_\d+\.tar$"),
]
# With --include-asdm: also fetch raw ASDMs (HUGE; only needed to recalibrate).
ASDM_PATTERNS = [
    re.compile(r"\.asdm\.sdm\.tar$"),
]
# With --include-products: also pull every other .pbcor.fits (calibrator imaging).
PRODUCT_PATTERNS = [
    re.compile(r"\.pbcor\.fits$"),
]


def filter_urls(urls, include_asdm=False, include_products=False,
                include_project_tar=False):
    keep = []
    for u_ in urls:
        if any(p.search(u_) for p in MS_PATTERNS):
            keep.append(u_)
        elif include_project_tar and any(p.search(u_) for p in PROJECT_TAR_PATTERNS):
            keep.append(u_)
        elif include_asdm and any(p.search(u_) for p in ASDM_PATTERNS):
            keep.append(u_)
        elif include_products and any(p.search(u_) for p in PRODUCT_PATTERNS):
            keep.append(u_)
    return keep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true",
                    help="Query archive and print what would be downloaded; do not fetch.")
    ap.add_argument("--proposals", nargs="*", default=None,
                    help="Only download these proposal_ids.")
    ap.add_argument("--sources", nargs="*", default=None,
                    help="Only download these source Names.")
    ap.add_argument("--match-radius", type=float, default=30.0,
                    help="Position-match radius in arcsec (default 30; mosaics may need 120).")
    ap.add_argument("--include-project-tar", action="store_true",
                    help="Fetch the full Member OUS product tarball "
                         "(_001_of_001.tar -- bundles every imaged spw and "
                         "every calibrator; can be hundreds of GB).")
    ap.add_argument("--include-asdm", action="store_true",
                    help="Also fetch the raw ASDM tarballs (HUGE: 50-500+ GB per MOUS). "
                         "Needed only if you want to re-run the calibration pipeline.")
    ap.add_argument("--max-gb-per-target", type=float, default=600.0,
                    help="Skip targets whose total download exceeds this size (default 600 GB).")
    ap.add_argument("--include-products", action="store_true",
                    help="Also fetch all .pbcor.fits files including calibrator imaging "
                         "(default already pulls *_sci.spw* science-target products).")
    ap.add_argument("--extract", action="store_true",
                    help="After download, extract each tarball.")
    ap.add_argument("--login-name", default=None,
                    help="ALMA archive username; only set for proprietary data.")
    ap.add_argument("--force", action="store_true",
                    help="Ignore the per-target .DONE marker and re-query / re-download.")
    args = ap.parse_args()

    # Lazy import so --help doesn't require astroquery
    from astroquery.alma import Alma

    # Load lists
    obs = Table.read(BEST_OBS_FITS)
    inv = Table.read(INVENTORY_CSV, format="csv")
    # FITS string columns come back as bytes (|S14); decode for comparisons.
    def _s(col):
        return np.array([
            x.decode() if isinstance(x, (bytes, bytearray)) else str(x)
            for x in col
        ])
    obs_pid  = _s(obs["proposal_id"])
    obs_name = _s(obs["Name"])
    on_disk  = {str(r["proposal_id"]) for r in inv if r["on_disk"] == "yes"}

    if args.force:
        mask = np.ones(len(obs), dtype=bool)
    else:
        mask = ~np.isin(obs_pid, list(on_disk))
    if args.proposals:
        mask &= np.isin(obs_pid, args.proposals)
    if args.sources:
        mask &= np.isin(obs_name, args.sources)
    todo = obs[mask]

    todo_pid  = obs_pid[mask]
    todo_name = obs_name[mask]
    print(f"Targets to download: {len(todo)} sources across "
          f"{len(set(todo_pid))} proposals")
    for i, r in enumerate(todo):
        pi = r['PI'].decode() if isinstance(r['PI'], (bytes, bytearray)) else r['PI']
        print(f"  {todo_pid[i]:16s}  {todo_name[i]:25s}  "
              f"PI={pi:25s}  res={float(r['best_res_arcsec']):.3f}\"")
    if args.dry_run and len(todo) == 0:
        return

    # Optional login (only needed for proprietary data)
    alma = Alma()
    if args.login_name and not args.dry_run:
        try:
            alma.login(args.login_name)
        except Exception as e:
            print(f"Login failed (continuing for public data only): {e}")

    UVDATA_ROOT.mkdir(parents=True, exist_ok=True)
    PROPOSAL_CACHE = {}

    for i, r in enumerate(todo):
        pid  = todo_pid[i]
        name = todo_name[i]
        # Find original source RA/Dec from full sample
        # (best_obs has dec; pull RA from sample table)
        ra_src, dec_src = lookup_radec(name)
        tgt_dir = UVDATA_ROOT / pid / name
        marker = tgt_dir / ".DONE"
        if marker.exists() and not args.force:
            print(f"[skip] {pid}/{name}  already DONE")
            continue
        if not args.dry_run:
            tgt_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n=== Querying {pid} / {name} (RA={ra_src:.4f}, Dec={dec_src:.4f}) ===")
        # 1) Pull all rows for the proposal once (cached per proposal).
        #    Alma.query(payload=...) returns 0 in newer astroquery versions;
        #    use TAP directly against ivoa.obscore.
        if pid not in PROPOSAL_CACHE:
            tap = alma.query_tap(
                "SELECT * FROM ivoa.obscore "
                f"WHERE proposal_id='{pid}'"
            )
            qp = tap.to_table() if hasattr(tap, "to_table") else tap
            PROPOSAL_CACHE[pid] = qp
            print(f"  proposal {pid}: {len(qp)} archive rows")
        qp = PROPOSAL_CACHE[pid]
        if len(qp) == 0:
            print(f"  proposal {pid} not found in archive; skipping")
            continue
        # 2) Filter to rows whose science target lies near our source
        qc = SkyCoord(np.asarray(qp["s_ra"], float)*u.deg,
                      np.asarray(qp["s_dec"], float)*u.deg)
        sep = qc.separation(
            SkyCoord(ra_src*u.deg, dec_src*u.deg)).arcsec
        q = qp[sep < args.match_radius]
        if len(q) == 0:
            print(f"  no archive rows within {args.match_radius}\" of {name}; "
                  f"skipping (proposal had {len(qp)} rows for other targets)")
            continue

        # Pick MOUS with smallest spatial_resolution = the row that delivered
        # the best resolution we identified earlier.
        if "spatial_resolution" in q.colnames:
            order = np.argsort(np.asarray(q["spatial_resolution"], float))
            q = q[order]
        mouses = []
        for row in q:
            mous = row["member_ous_uid"] if "member_ous_uid" in q.colnames else row["Member ous id"]
            if mous not in mouses:
                mouses.append(mous)

        urls_to_get = []
        for mous in mouses:
            try:
                info = alma.get_data_info(mous, expand_tarfiles=True)
            except Exception as e:
                print(f"  get_data_info({mous}) failed: {e}")
                continue
            urls = list(info["access_url"])
            urls = filter_urls(urls,
                               include_asdm=args.include_asdm,
                               include_products=args.include_products,
                               include_project_tar=args.include_project_tar)
            print(f"  MOUS {mous}: {len(urls)} URLs after filtering")
            urls_to_get.extend(urls)

        urls_to_get = sorted(set(urls_to_get))
        if not urls_to_get:
            print("  no MS-bearing URLs matched; nothing to fetch")
            continue

        # Estimate total size from content_length on the data_info rows we cached
        size_b = 0
        for mous in mouses:
            try:
                info = alma.get_data_info(mous, expand_tarfiles=True)
                for u_, s in zip(info["access_url"], info["content_length"]):
                    if u_ in urls_to_get and s:
                        size_b += int(s)
            except Exception:
                pass
        gb = size_b / 1e9
        print(f"  {len(urls_to_get)} URLs, total {gb:.2f} GB")
        if gb > args.max_gb_per_target and not args.dry_run:
            print(f"  REFUSING to auto-download {gb:.0f} GB "
                  f"(>{args.max_gb_per_target} GB). Re-run with "
                  f"--max-gb-per-target {int(gb)+1} to bypass.")
            continue

        if args.dry_run:
            for u_ in urls_to_get[:5]:
                print(f"    would fetch: {u_}")
            if len(urls_to_get) > 5:
                print(f"    ... ({len(urls_to_get)-5} more)")
            continue

        alma.cache_location = str(tgt_dir)
        try:
            local = alma.download_files(urls_to_get,
                                         savedir=str(tgt_dir),
                                         continuation=True)
            print(f"  fetched {len(local)} files into {tgt_dir}")
        except Exception as e:
            print(f"  download failed: {e}")
            continue

        if args.extract:
            ms_found = extract_tarballs(tgt_dir, [str(p) for p in local])
            if not ms_found:
                print(f"  No calibrated MS yet -- run scriptForPI.py to "
                      f"generate one from the raw ASDM:")
                # find the scriptForPI in the extracted tree
                import glob
                for s in glob.glob(str(tgt_dir / "**" / "scriptForPI.py"),
                                   recursive=True):
                    print(f"    {s}")
        marker.write_text("DONE\n")


def _bcol(col):
    out = []
    for x in col:
        if isinstance(x, (bytes, bytearray)):
            x = x.decode()
        out.append(str(x).strip())
    return np.array(out)


def extract_tarballs(tgt_dir, tar_paths):
    """Extract every *.tar in tar_paths into tgt_dir; return True if any
    .ms / .ms.split.cal directory is created."""
    import subprocess
    found_ms = False
    for tp in tar_paths:
        if not tp.endswith(".tar"):
            continue
        if not os.path.isfile(tp):
            continue
        print(f"    extracting {os.path.basename(tp)}")
        try:
            subprocess.run(["tar", "xf", tp, "-C", str(tgt_dir)],
                           check=True)
        except subprocess.CalledProcessError as e:
            print(f"      extract failed: {e}")
            continue
    # search for any MS directory
    import glob
    for pat in ("**/*.ms", "**/*.ms.split.cal", "**/*calibrated*.ms"):
        if glob.glob(os.path.join(str(tgt_dir), pat), recursive=True):
            found_ms = True
            break
    return found_ms


def lookup_radec(name):
    for path in ("data/myso_sample_2026_full.fits",
                 "data/myso_sample_2026_kept.fits",
                 "data/myso_sample_2026_deep.fits"):
        if not os.path.exists(path):
            continue
        full = Table.read(path)
        names = _bcol(full["Name"])
        m = names == name
        if m.any():
            r = full[m][0]
            return float(r["RAJ2000"]), float(r["DEJ2000"])
    raise KeyError(f"{name} not in any catalog")


if __name__ == "__main__":
    main()
