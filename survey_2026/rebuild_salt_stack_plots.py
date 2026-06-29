"""Retroactively rebuild per-source nacl/kcl/joint stack PNGs from existing
.spec.npz files. Re-derives the 35Cl v=0/1/2 line list per source by
listing matching .spec.npz files in source_<NN>/, re-stacks them with the
v_LSR from the per-target catalog, and writes:

  nacl_stack.png
  kcl_stack.png
  naclkcl_stack.png
  naclkcl_combined_stack.png (+ .npz)

with the labels showing which lines went into each stack.

Run this once after the line_pipeline upgrade to refresh older outputs.
"""
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

import line_pipeline as lp

ROOT = Path("/orange/adamginsburg/salt/survey_2026")
ANALYSIS = ROOT / "analysis_products"

ISOTOP_RE = re.compile(r"^(NaCl|KCl)_v[012]_")


def vlsr_for(target):
    if target in lp.__dict__.get("MANUAL_VLSR", {}):
        return lp.MANUAL_VLSR[target]
    # Fall back to data/rms.fits
    rms = ROOT / "data/rms.fits"
    if rms.exists():
        from astropy.table import Table
        t = Table.read(rms)
        names = np.array([str(n) if not isinstance(n, bytes) else n.decode()
                          for n in t["Name"]])
        m = np.where(names == target)[0]
        if len(m):
            v = float(t["vLSR"][m[0]])
            if np.isfinite(v):
                return v
    # Hand overrides for non-RMS targets
    manual = {"MonR2-IRS3": 10.0, "MonR2-IRS2": 10.0,
              "NGC6334I": -7.0, "NGC6334IN": -3.5,
              "Orion_SrcI": 5.0, "Orion-SrcI": 5.0}
    return manual.get(target, None)


def per_source_stack(src_dir, vlsr_kms):
    """Re-derive {species: dict(vaxis, spec, sigma, n_lines, lines_used)}
    from the .spec.npz files in src_dir."""
    out = {}
    for species in ("NaCl", "KCl"):
        candidates = sorted([p for p in src_dir.glob(f"{species}_v[012]_*.spec.npz")
                              if ISOTOP_RE.match(p.stem.replace(".spec", ""))])
        if len(candidates) < 2:
            continue
        v_grid = np.arange(-lp.VWIN_KMS, lp.VWIN_KMS + 0.5, 1.0)
        stacks, sigmas, used = [], [], []
        for npz_path in candidates:
            arr = np.load(npz_path)
            if "vaxis" not in arr.files or "spec" not in arr.files:
                continue
            v_shift = arr["vaxis"] - vlsr_kms
            order = np.argsort(v_shift)
            interp = np.interp(v_grid, v_shift[order], arr["spec"][order],
                                left=np.nan, right=np.nan)
            if not np.isfinite(interp).any():
                continue
            stacks.append(interp)
            sigmas.append(float(arr["sigma"]))
            used.append(npz_path.stem.replace(".spec", ""))
        if len(stacks) < 2:
            continue
        weights = 1.0 / np.array(sigmas) ** 2
        stack = np.nansum(np.array(stacks) * weights[:, None], axis=0) / np.nansum(weights)
        sigma = 1.0 / np.sqrt(np.nansum(weights))
        if not np.isfinite(stack).any() or np.nanstd(stack) == 0:
            continue
        out[species] = dict(vaxis=v_grid, spec=stack, sigma=sigma,
                              n_lines=len(stacks), lines_used=used)
    return out


def main():
    n_done = 0
    for src_dir in sorted(ANALYSIS.glob("*/2*/source_*")):
        target = src_dir.parent.parent.name
        proposal = src_dir.parent.name
        m = re.match(r"^source_(\d{2})$", src_dir.name)
        if not m:
            continue
        sid = int(m.group(1))
        vlsr = vlsr_for(target)
        if vlsr is None:
            print(f"  skip {target}/{proposal}/src{sid:02d}: no vlsr")
            continue
        stacks = per_source_stack(src_dir, vlsr)
        if not stacks:
            continue
        lp.plot_salt_stack(stacks, src_dir, sid)
        # Also overwrite the legacy NaCl_stack.npz / KCl_stack.npz with the
        # lines_used field
        for sp, d in stacks.items():
            np.savez(src_dir / f"{sp}_stack.npz",
                      vaxis=d["vaxis"], spec=d["spec"], sigma=d["sigma"],
                      n_lines=d["n_lines"],
                      lines_used=np.array(d["lines_used"]))
        n_done += 1
        species_str = "/".join(sorted(stacks.keys()))
        print(f"  {target}/{proposal}/src{sid:02d}: {species_str} "
              f"(v_LSR={vlsr:+.1f})")
    print(f"\nupdated {n_done} sources")


if __name__ == "__main__":
    main()
