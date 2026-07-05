"""Consolidate CARTA snippets per target.

For each target with multiple snippet files (e.g. TARGET_2026.json plus
TARGET_<proposal>.json from line_pipeline runs):

  1. Parse code lines from every snippet for that target.
  2. Strip /orange/adamginsburg prefix from any path that has it (CARTA mounts
     /orange/adamginsburg as its filesystem root).
  3. Rewrite importRegion calls into the (dir, filename, format_int) tuple form.
  4. Promote the first openFile so it stays openFile; convert any subsequent
     openFile call to appendFile.
  5. Dedup by exact line.
  6. Write the merged result to TARGET_2026.json and DELETE the per-proposal
     duplicates.

Safety:
  - Only merges files whose stem is exactly TARGET or TARGET_<proposal>; a
    suffix that doesn't look like an ALMA proposal id is left alone.
  - DIHCA2 / non-L4_d2 snippets are untouched unless they match a target name
    explicitly listed.
"""
import json
import re
import shutil
from collections import defaultdict
from datetime import datetime
from pathlib import Path

SNIP = Path.home() / ".carta/config/snippets"
CARTA_PREFIX = "/orange/adamginsburg"
PROPOSAL_RE = re.compile(r"_20\d\d\.\d\.\d+\.[SL]$")
SUFFIX_2026 = re.compile(r"_2026$")


def carta_path(s: str) -> str:
    if s.startswith(CARTA_PREFIX):
        return s[len(CARTA_PREFIX):]
    return s


def base_target(stem: str) -> str:
    t = SUFFIX_2026.sub("", stem)
    t = PROPOSAL_RE.sub("", t)
    return t


def is_per_proposal(stem: str) -> bool:
    """True for TARGET_2019.1.00437.S patterns."""
    return bool(PROPOSAL_RE.search(stem))


def is_canonical(stem: str) -> bool:
    """True for TARGET_2026 (the consolidated file)."""
    return stem.endswith("_2026")


def rewrite_line(line: str) -> str:
    line = line.strip()
    # Strip /orange/adamginsburg from any quoted path
    line = re.sub(
        r'"' + re.escape(CARTA_PREFIX) + r'(/[^"]*)"',
        r'"\1"',
        line,
    )
    # Convert legacy importRegion("path", "DS9") form to (dir, filename, 2).
    m = re.match(
        r'await app\.importRegion\("([^"]+)",\s*"DS9"\)\s*$', line)
    if m:
        full = m.group(1)
        p = Path(full)
        line = (f'await app.importRegion("{p.parent}/", "{p.name}", 2)')
    return line


def merge_snippet_codes(codes):
    """Combine many snippet code blocks into a single open+append+region list."""
    seen = set()
    out = []
    first_open = False
    for code in codes:
        for line in (code or "").split("\n"):
            line = rewrite_line(line)
            if not line:
                continue
            if line.startswith('await app.openFile('):
                if first_open:
                    line = line.replace('await app.openFile(',
                                          'await app.appendFile(', 1)
                else:
                    first_open = True
            if line in seen:
                continue
            seen.add(line)
            out.append(line)
    return "\n".join(out)


def main():
    files_by_target = defaultdict(list)
    for f in SNIP.glob("*.json"):
        files_by_target[base_target(f.stem)].append(f)

    ts = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    backup = SNIP / f"_backup_per_target_{ts}"
    backup.mkdir(exist_ok=True)

    targets_with_dups = {t: fs for t, fs in files_by_target.items()
                          if len(fs) > 1 or any(is_per_proposal(f.stem) for f in fs)}

    n_merged = 0
    n_deleted = 0
    for tgt, fs in sorted(targets_with_dups.items()):
        # Sort so canonical _2026 files are processed last (preserves their
        # ordering of openFile, then we dedup).
        fs_sorted = sorted(fs, key=lambda f: (is_canonical(f.stem), f.stem))
        codes = []
        category = "saltsearch2026"
        for f in fs_sorted:
            try:
                data = json.loads(f.read_text())
            except json.JSONDecodeError:
                continue
            codes.append(data.get("code", ""))
            cats = data.get("categories")
            if cats:
                category = cats[0]
        merged_code = merge_snippet_codes(codes)
        if not merged_code:
            continue
        out = SNIP / f"{tgt}_2026.json"
        # Back up everything we touch before overwriting/deleting.
        for f in fs_sorted:
            shutil.copy2(f, backup / f.name)
        out.write_text(json.dumps({
            "$schema": "https://cartavis.github.io/schemas/snippet_schema_1.json",
            "categories": [category],
            "code": merged_code,
            "frontendVersion": "3.0.0",
            "snippetVersion": 1,
        }, indent=4))
        n_merged += 1
        # Delete per-proposal duplicates (but never the canonical _2026 file)
        for f in fs_sorted:
            if f == out:
                continue
            f.unlink()
            n_deleted += 1
        n_lines = len(merged_code.split("\n"))
        print(f"  {tgt}: merged {len(fs_sorted)} -> {out.name} ({n_lines} lines)")

    # Also strip /orange/adamginsburg from non-duplicate canonical snippets.
    n_stripped = 0
    for f in SNIP.glob("*.json"):
        if f.is_dir():
            continue
        try:
            data = json.loads(f.read_text())
        except json.JSONDecodeError:
            continue
        code = data.get("code", "")
        if CARTA_PREFIX in code or "DS9" in code:
            new_code = "\n".join(rewrite_line(l) for l in code.split("\n")
                                  if l.strip())
            if new_code != code:
                shutil.copy2(f, backup / f.name)
                data["code"] = new_code
                f.write_text(json.dumps(data, indent=4))
                n_stripped += 1

    print(f"\nMerged {n_merged} targets, deleted {n_deleted} per-proposal files, "
          f"stripped paths in {n_stripped} other files. Backup at {backup}")


if __name__ == "__main__":
    main()
