"""
Inventory which proposal codes (from the salt MYSO target list) already have
pipeline data on disk somewhere under /orange/adamginsburg.

Output: data/myso_alma_disk_inventory.csv  (and prints a summary).

Each row in myso_alma_best_obs.fits has a proposal_id (the project that
delivered the highest-resolution observation of that source).  We resolve
each proposal_id to a directory tree on disk if one exists, by:

  1. Globbing /orange/adamginsburg/salt/**/<proposal_id>/  (depth<=4)
  2. Globbing /orange/adamginsburg/**/<proposal_id>/       (depth<=3)
  3. Falling back to a manual override map for projects we know we have
     data for under non-standard paths (e.g. miriam/W33A/...).

Whatever isn't found is what download_alma_data.py needs to fetch.
"""

import os
import glob
from pathlib import Path
from astropy.table import Table


# Known on-disk locations for proposals we already have but stored under
# non-standard paths (proposal_id is not in the directory name).
MANUAL_OVERRIDE = {
    # AFGL2591 / MAUDS / Maud+ — already in W33A reduction tree
    "2018.1.00458.S": [
        "/orange/adamginsburg/salt/miriam/W33A/W33_mauds2/2018.1.00458.S",
        "/orange/adamginsburg/salt/DIHCA_salts/W33_mauds2/2018.1.00458.S",
    ],
    "2022.1.00700.S": [
        "/orange/adamginsburg/salt/sanhueza/2022.1.00700.S",
    ],
}

# Restricted to /orange/adamginsburg/salt only -- this is where ALMA project
# data lives.  /orange/adamginsburg and /blue/adamginsburg roots contain too
# many irrelevant subdirs (conda envs, ALMA_IMF, etc.) to scan exhaustively.
SEARCH_ROOTS = [
    ("/orange/adamginsburg/salt", 6),
]
# Skip these subtrees -- nothing useful + slow to walk
PRUNE_DIRS = (
    "archive", "code", "notebooks", "papers", "scratch",
    ".ipynb_checkpoints", ".git", ".cache", "survey_2026",
    "dihca2_source_data", "dihca2", "dihca",  # different sample
)
_DIR_CACHE = None


def _scan_all_dirs():
    """One-shot find: list every directory whose basename looks like an
    ALMA project code (e.g. 2022.1.01344.S or 2019.1.00195.L).
    Cached so repeated calls are free."""
    global _DIR_CACHE
    if _DIR_CACHE is not None:
        return _DIR_CACHE
    cands = []
    for root, depth in SEARCH_ROOTS:
        if not os.path.isdir(root):
            continue
        prune = " -o ".join(f"-name {d}" for d in PRUNE_DIRS)
        cmd = (
            f'find {root} -maxdepth {depth} '
            f'\\( \\( {prune} \\) -prune \\) -o '
            f'-type d -regextype posix-extended '
            f'-regex ".*/[0-9]{{4}}\\.[0-9]\\.[0-9]+\\.[SL]$" -print '
            f'2>/dev/null'
        )
        out = os.popen(cmd).read()
        cands.extend(l.strip() for l in out.splitlines() if l.strip())
    print(f"  (scanned {len(cands)} candidate proposal dirs)", flush=True)
    _DIR_CACHE = cands
    return cands


def find_on_disk(pid):
    if pid in MANUAL_OVERRIDE:
        hits = [p for p in MANUAL_OVERRIDE[pid] if os.path.isdir(p)]
        if hits:
            return hits
    return sorted(p for p in _scan_all_dirs() if os.path.basename(p) == pid)


def main():
    t = Table.read("data/myso_alma_best_obs.fits")
    pids = sorted(set(t["proposal_id"]))
    rows = []
    for i, pid in enumerate(pids, 1):
        sources = [r["Name"] for r in t if r["proposal_id"] == pid]
        print(f"[{i:2d}/{len(pids)}] {pid} ...", flush=True)
        hits = find_on_disk(pid)
        print(f"           -> {'ON_DISK' if hits else 'no'}  {hits}", flush=True)
        rows.append({
            "proposal_id": pid,
            "n_sources": len(sources),
            "sources": ";".join(sources),
            "on_disk": "yes" if hits else "no",
            "paths": ";".join(hits) if hits else "",
        })
    out = Table(rows=rows)
    out.write("data/myso_alma_disk_inventory.csv", format="csv", overwrite=True)
    on = sum(1 for r in rows if r["on_disk"] == "yes")
    print(f"On disk      : {on}/{len(rows)} proposals")
    print(f"To download  : {len(rows)-on}/{len(rows)}")
    print()
    print(f"{'proposal':16s} {'n':>3s} {'src':6s} paths")
    for r in rows:
        print(f"  {r['proposal_id']:14s} {r['n_sources']:3d} {r['on_disk']:6s} {r['paths']}")


if __name__ == "__main__":
    main()
