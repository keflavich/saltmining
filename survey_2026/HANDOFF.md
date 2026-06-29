# survey_2026 handoff doc

Date: 2026-06-29. Paper at `/orange/adamginsburg/salt/demography_2026/` (Overleaf git-synced). Code at `/orange/adamginsburg/salt/survey_2026/`.

## Goal

Multi-month ALMA salt-search demography paper: NaCl/KCl/H2O detections+ULs across MYSO disks at d ≤ 2 kpc, L ≥ 10⁴ L⊙. 48-source sample drawn from RMS + uchii_2026 cross-match.

## Paper status

- `main.pdf` compiles (660KB) under pdflatex (no engine other than that — load `module load texlive`; use `/apps/texlive/2023/bin/x86_64-linux/pdflatex`).
- Tables converted to portrait via `\startlongtable` (NO `longrotatetable`). User insisted strongly.
- Citation references inside tables go to numbered `\tablerefs{[N] \citet{X}; ...}` footnotes — not inline. See `build_target_table.py:citet_to_num` for the pattern.
- Figures: 10 paper figures inserted as `\input{}`-equivalent `\begin{figure*}` blocks in `main.tex § Figures`.

## Critical files

```
survey_2026/
├── data/
│   ├── sources_L4_d2.csv           # 48-row master sample (any add: also update build_target_table.py)
│   ├── literature_detections.csv   # per-target lit overrides (Ginsburg+19, Plambeck+13, ...)
│   ├── vlsr_from_literature.json
│   ├── vlsr_from_data.json         # vlsr measured from corroborating lines
│   ├── per_target_paper.csv        # auto-built per-target detection table
│   ├── evidence_audit.csv          # per-source warning flags
│   ├── detected_lines_report.{csv,md}      # snr>=3.5 detected lines per source
│   └── detected_lines_with_contaminants.{csv,md}  # XCLASS-flagged contaminants
├── analysis_products/<target>/<proposal>/
│   ├── continuum_sources.csv       # mm cont sources (id, ra, dec, peak_Jybeam, snr)
│   ├── line_measurements.csv       # per (source, line) snr/peak/integ/peak_v
│   ├── source_<NN>/                # NaCl_stack.npz, kcl_stack, nacl/kcl/joint stack pngs
│   │   ├── *_diagnostic.png        # per-line diagnostic
│   │   ├── source_<NN>_<line>_mom0.fits  # mom0
│   │   ├── source_<NN>_<line>_mom1.fits
│   │   └── *.spec.npz              # extracted spectrum (vaxis, spec, sigma)
│   ├── spectrum_panels_src<NN>.png # N-panel SPW spectrum w/ lineid labels
│   ├── kinematic_stack/aligned_by_<guide>.{png,npz}
│   └── by_type/                    # symlink tree grouping products by class
├── lineid/                         # XCLASS line-id scaffolding
│   ├── <target_label>/
│   │   ├── *_band6.molfit          # CH3OH/CH3OCHO/CH3CN/NaCl/KCl/H2O active by default
│   │   ├── spectra/<label>_spw*.fits  # source-pixel spectra in K
│   │   ├── models/spw*.model.npz   # XCLASS LTE model
│   │   └── plots/spw*.png          # data + model overlay
│   ├── extract_source_spectra.py
│   ├── xclass_runner.py
│   └── run_xclass.sh               # activates the xclass env + LD_LIBRARY_PATH
├── build_*.py                      # per-table builders (target/multisurvey/per_target/...)
├── line_pipeline.py                # source detection + line measurement entry point
├── run_line_pipeline_slurm.sh      # slurm wrapper
├── launch_l4_d2_pipeline.py        # auto-picks best proposal per target + submits
├── make_paper_figures.py
├── HANDOFF.md (this file)
└── memory/                         # auto-memory (claude-private; ~/.claude/projects/...)
```

## Pipelines

### Line measurement: `line_pipeline.py`
- Inputs: `--proposal <id> --target <name> --vlsr <kms> --distance-kpc <d> [--on-kms <kms>]`
- Outputs `analysis_products/<target>/<proposal>/`. Writes continuum_sources.csv, line_measurements.csv, per-source diagnostic PNGs, NaCl/KCl stacks.
- Submit via `sbatch run_line_pipeline_slurm.sh --proposal ...`. QOS rule: heavy jobs use `--account=astronomy-dept --qos=astronomy-dept-b` (the slurm script sets this).
- Bulk launch: `python launch_l4_d2_pipeline.py [--force]` iterates all L4_d2 sources, picks best downloaded proposal, submits.

### XCLASS line ID: `lineid/run_xclass.sh <label>`
- Activates `/blue/adamginsburg/adamginsburg/miniconda3/envs/xclass`
- Sets `LD_LIBRARY_PATH=/apps/libgfortran/3.0.0/lib:$LD_LIBRARY_PATH` (needed by `myNewXCLASS.exe`)
- Runs `xclass_runner.py <label>` which calls `task_myXCLASS.myXCLASS()` per SPW
- Tested working for G326.6618+00.5207, G015.0357, G345.5043
- Output: `lineid/<label>/{models, plots}/spw*.{npz,png}`
- Transition table from XCLASS run is at `/orange/adamginsburg/software/XCLASS-Interface/run/myXCLASS/job__<ts>__<rand>/transition_energies.dat`

### Contaminant ID: `find_contaminants_from_xclass.py`
- For each suspect 5σ detection, search the XCLASS transition_energies.dat for transitions whose Doppler-shifted freq is within ±25 MHz of the observed peak freq
- Output: `data/detected_lines_with_contaminants.{csv,md}`
- Key finding so far: G015.0357 NaCl_v2_J19-18 7.4σ likely CH3OCHO; G326 H2O_v2 candidate CH3OCHO blend

## Recent definitive paper findings

### Confirmed salt detections
- **Orion-SrcI** (Ginsburg+2019): NaCl/KCl/H2O/SiO/SO + COMs incontrovertible (literature override)
- **MonR2-IRS3** (G213.71-12.60): NaCl 4-line stack snr=4.8σ peak=2.07 mK at src01 (2019.1.00437.S)

### Strongly suspect "detections" needing XCLASS vetting
- **G015.0357 NaCl_v2_J19-18 7.4σ** (M17 SW UCHII): v=2 detection without v=0 corroboration; CH3OCHO blend likely (XCLASS predicts 10.7 K at 243.5646 GHz, 5.5 MHz from peak)
- **G326.6618 H2O_v2 16.3σ + H2O 5.3σ**: both H2O peaks at v_obs = -49 km/s (10 km/s blue of vsys = -39.6); CH3OCHO blend predicted
- **G345.5043 NaCl_v2_J17-16 6.6σ**: same v=2-only pattern; needs XCLASS check
- **MonR2-IRS3 H30α detection in 2022.1.00285.S**: peak_v=+42 km/s (32 km/s offset!) — almost certainly not RRL

### Confirmed non-detections (tight ULs)
- 8 sources have TIGHT ULs at <3 mK NaCl, <5 mK H2O: NGC6334IN, G081.8789, G081.6802B, G345.4938, G342.1156, G345.0061C, G189.0307, MonR2-IRS3 (paradox — needs handling)

## Open tasks (62 total)

Use `TaskList` to view all. Major categories:

| Category | Count | Notes |
|---|---|---|
| Self-consistency vs lit | 8 | #26-33 |
| Download + analyze new sources | 16 | #34-49, #56-57, #59 |
| ALMA-blocked (dec>+60) | 7 | #50-55, #58 — close as not-doable |
| Partial analysis cleanup | 3 | #60-62 |

### Most urgent
- **#27 NGC6334I-mm1b source ID**: brightest src05 in pipeline ≠ literature mm1b. 7 species disagreements with detections.tex.
- **#33 Auto-build detections.tex** from per_target + lit (root cause fix — eliminates drift).
- **#35 OrionB-Flame re-imaging**: 2017.1.01102.S only delivered mfs continuum; cubes need CASA tclean from raw vis.
- **In-flight slurm batch** (jobs 35934221-35934242): 22 line_pipeline jobs PD on QOSGrpCpuLimit — will run as slots free. NGC6334I (12 proposals), NGC6334IN (4), G268.4222 (80 cubes!), G232.6207 (2), G010.8411 (4).

### Slurm batch monitoring
```bash
squeue -u adamginsburg --name=line_pipe -h | wc -l   # count
squeue -u adamginsburg --name=line_pipe -h -t R | wc -l  # running
ls logs/line_pipe_*.log | head -5      # recent log files
grep -l Traceback logs/line_pipe_3593*.log    # find failures
```

After batch completes:
1. `python build_evidence_audit.py`
2. `python build_per_target_paper_table.py`
3. `python build_data_summary.py`
4. `python build_byproduct_symlinks.py`
5. `python make_paper_figures.py`
6. Compile main.pdf

## Operational rules (user-established)

1. **Caveman mode active**: drop articles/filler/pleasantries. Fragments OK.
2. **No mock databases in tests**: use real cubes.
3. **No bare `except`** in Python — only specific exception types per CLAUDE.md.
4. **Don't use 'honest'**.
5. **CARTA paths NOT include /orange/adamginsburg**: CARTA mounts that as root. Use `/salt/...` paths. importRegion syntax: `(dir, filename, 2)` (3-arg tuple).
6. **SLURM heavy jobs**: `--account=astronomy-dept --qos=astronomy-dept-b`. Never astronomy-dept QOS alone. Light jobs: adamginsburg/adamginsburg.
7. **No bulk deletes** (no `rm -rf`, no wildcard deletes). Overwrite in place; backups for renames.
8. **Don't cancel running slurm jobs** unless asked.
9. **TMPDIR=/blue/adamginsburg/adamginsburg/tmp** for spectral_cube intermediates (/tmp fills up).
10. **All tables PORTRAIT** (no `\rotate`, no landscape). User strong preference.

## Recent commits (latest first)
```
443bf15  Doc OrionB-Flame DL: only mfs continuum delivered
2422182  OrionB-Flame DL: corrected coords + 502 retry loop
f6ae438  Queue OrionB-Flame 2017.1.01102.S download
603b822  Builders aligned with LaTeX fixes
9967951  Fix stray nodata + obs_params single-unit + upper_limits output-loop
04b9677  v_LSR column added to targets table + figures section
b5bee02  Portrait-mode tables + targets meta + ALMA companion split
2917e95  by_type/ symlink tree
fa3da6e  Multi-survey cross-match table
d50b76e  CARTA snippet WCS-filter + stack FITS files
```

## Key invariants

- Every `analysis_products/<target>/<proposal>/source_<NN>/*_diagnostic.png` should be named `source_<NN>_<line>_diagnostic.png` (bulk-renamed 1676 files for this).
- `data/sources_L4_d2.csv` is the canonical 48-target sample. Adding sources requires also updating PARENT_REGION in `build_multisurvey_table.py` and SPECIAL_DIST_REF / SPECIAL_LUM_REF in `build_target_table.py`.
- `literature_detections.csv` takes precedence over pipeline-derived per_target rows.
- CARTA snippets at `~/.carta/config/snippets/<target>_2026.json` use `/salt/...` paths and 3-arg `importRegion`.

## Failure modes seen

1. **ALMA archive 502** on `Alma.get_data_info` — retry loop in download_orionb_flame.py
2. **Splatalogue 503** — use XCLASS transition table instead (`find_contaminants_from_xclass.py`)
3. **myNewXCLASS.exe `libgfortran.so.3 not found`** — set `LD_LIBRARY_PATH=/apps/libgfortran/3.0.0/lib`
4. **deluxetable+startlongtable "Output loop 200 dead cycles"** — table too wide; reduce columns or font, do NOT keep `\startlongtable` if table is too narrow
5. **r'\\nodata' rendering as line-break + "nodata" text**: use `r'\nodata'` (one backslash)
6. **find_misplaced_cubes**: many uvdata/<proposal>/<target>/ dirs contain off-target cubes (multi-pointing projects). FOV filter in build_snippets_recent excludes them from CARTA but uvdata not pruned.

## Where things are stored

- Paper repo: `/orange/adamginsburg/salt/demography_2026/` (git-tracked, Overleaf)
- Code repo: `/orange/adamginsburg/salt/survey_2026/` (git-tracked)
- ALMA cubes: `/orange/adamginsburg/salt/survey_2026/uvdata/<proposal>/<target>/`
- Pipeline outputs: `/orange/adamginsburg/salt/survey_2026/analysis_products/`
- Slurm logs: `/orange/adamginsburg/salt/survey_2026/logs/line_pipe_<jobid>.log`
- Memory (claude private): `~/.claude/projects/-orange-adamginsburg-salt/memory/`
- XCLASS job artifacts: `/orange/adamginsburg/software/XCLASS-Interface/run/myXCLASS/job__<ts>__<rand>/`
- Lit ref source: `/orange/adamginsburg/salt/Orion_ALMA_2016.1.00165.S/analysis/lines.py` (`disk_lines` + `absorbers` = 161 transitions)
- W51 XCLASS reference: `/orange/adamginsburg/w51/hotcores/linecatalog/ethanbhula/`

## Compaction note

This session ended after 60+ task creations and many edits. Memory entries pre-existing:
- `survey_2026 marked COMPLETE` was wrong — still active. User has been adding tasks (XCLASS, literature reconciliation, etc.).
- See `~/.claude/projects/-orange-adamginsburg-salt/memory/MEMORY.md` for the rest.
