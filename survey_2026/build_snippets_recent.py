"""Build CARTA snippets for recently-downloaded regions.

For each target, find local cont + cube files across all downloaded proposals
and write a snippet to ~/.carta/config/snippets/ under category saltsearch2026.
"""
import json
import re
from pathlib import Path

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
UVDIR = ROOT / "uvdata"
SNIP = Path.home() / ".carta/config/snippets"
SNIP.mkdir(parents=True, exist_ok=True)

CATEGORY = "saltsearch2026"

TARGETS = None  # set None to iterate every L4_d2 target with downloaded data


def normalize(s):
    return re.sub(r"[\s_\-+]", "", str(s)).lower()


def name_variants(s):
    base = str(s)
    out = {normalize(base),
           normalize(base.replace("+", "p").replace("-", "m")),
           normalize(base.replace("+", "").replace("-", ""))}
    return {v for v in out if v}


def find_files(target):
    """Return (best_cont, sorted_cubes) across all proposals for this target."""
    tgt_keys = name_variants(target)
    cont_paths, cube_paths = [], []
    for prop_dir in sorted(UVDIR.iterdir()):
        if not prop_dir.is_dir():
            continue
        for sub in prop_dir.iterdir():
            if not sub.is_dir():
                continue
            if normalize(sub.name) not in tgt_keys:
                continue
            for f in sub.iterdir():
                if not f.name.endswith(".I.pbcor.fits"):
                    continue
                if ".cont." in f.name:
                    cont_paths.append(f)
                elif ".cube." in f.name:
                    cube_paths.append(f)
    # best cont: multi-spw combined preferred
    cont_best = None
    if cont_paths:
        multi = [p for p in cont_paths if re.search(r"spw\d+(_\d+){2,}", p.name)]
        pool = multi if multi else cont_paths
        cont_best = sorted(pool, key=lambda p: p.stat().st_size)[0]
    return cont_best, sorted(cube_paths, key=lambda p: p.name)


CARTA_ROOT_PREFIX = "/orange/adamginsburg"


def carta_path(p):
    s = str(p)
    return s[len(CARTA_ROOT_PREFIX):] if s.startswith(CARTA_ROOT_PREFIX) else s


def find_region_files(target):
    """Return sorted list of *.reg files under analysis_products/<target>/<prop>/."""
    base = ROOT / "analysis_products" / target
    if not base.is_dir():
        return []
    regs = []
    for prop_dir in sorted(base.glob("2*")):
        for r in prop_dir.glob("*.reg"):
            regs.append(r)
    return regs


def build_snippet_code(cont, cubes, regs):
    lines = []
    rest = list(cubes)
    if cont is not None:
        lines.append(f'await app.openFile("{carta_path(cont)}")')
    elif rest:
        first = rest.pop(0)
        lines.append(f'await app.openFile("{carta_path(first)}")')
    for c in rest:
        lines.append(f'await app.appendFile("{carta_path(c)}")')
    for r in regs:
        rdir = carta_path(r.parent) + "/"
        rname = r.name
        lines.append(f'await app.importRegion("{rdir}", "{rname}", 2)')
    return "\n".join(lines)


def main():
    if TARGETS is None:
        import pandas as pd
        df = pd.read_csv(ROOT / "data/sources_L4_d2.csv")
        targets = sorted(set(df["name"]))
    else:
        targets = TARGETS
    for tgt in targets:
        cont, cubes = find_files(tgt)
        regs = find_region_files(tgt)
        if cont is None and not cubes:
            print(f"  {tgt}: NO LOCAL DATA — skipping")
            continue
        code = build_snippet_code(cont, cubes, regs)
        snippet = {
            "$schema": "https://cartavis.github.io/schemas/snippet_schema_1.json",
            "categories": [CATEGORY],
            "code": code,
            "frontendVersion": "3.0.0",
            "snippetVersion": 1,
        }
        out = SNIP / f"{tgt}_2026.json"
        out.write_text(json.dumps(snippet, indent=4))
        print(f"  wrote {out.name}: cont={'yes' if cont else 'no'}, "
              f"cubes={len(cubes)}, regs={len(regs)}")


if __name__ == "__main__":
    main()
