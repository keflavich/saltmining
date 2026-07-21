# DIHCA2 Imaging Plan - Quick Start Guide

## Overview

The imaging plan has been successfully generated! Here's what was created and how to use it.

## Generated Files

### 1. **dihca2_imaging_plan.json** (815 KB)
The main output file containing 844 CASA tclean imaging commands for 211 source groups across 65 fields.

### 2. **imaging_plan.py** (9.6 KB)
The main script that generates the imaging plan from the source catalog.

**Key features:**
- Reads `dihca_source_catalog.fits` (533 sources)
- Groups sources spatially using hierarchical clustering (5.12" threshold)
- Adaptively sizes images (512×512 or 1024×1024 pixels)
- Generates commands for all 4 spectral windows

**Regenerate the plan:**
```bash
python imaging_plan.py
```

### 3. **query_imaging_plan.py** (7.6 KB)
Utility script to query and analyze the imaging plan.

**Examples:**
```bash
# Show overall statistics
python query_imaging_plan.py --stats

# Find all groups for a specific field
python query_imaging_plan.py --field G008.672-00.682 --verbose

# Find groups with many sources
python query_imaging_plan.py --min-sources 10

# Find all 1024×1024 pixel images
python query_imaging_plan.py --imsize 1024

# Export commands for a specific field
python query_imaging_plan.py --field G345.493+01.469 --export
```

### 4. **run_imaging.py** (8.7 KB)
Script to execute the imaging commands in CASA.

**Usage in CASA:**
```bash
# Dry run to see what would be executed
casa --nologger -c "execfile('run_imaging.py'); run_field('G008.672-00.682', dry_run=True)"

# Actually run imaging for a specific field
casa --nologger -c "execfile('run_imaging.py'); run_field('G008.672-00.682')"

# Run a specific group
casa --nologger -c "execfile('run_imaging.py'); run_group('G008.672-00.682_group1')"
```

Or use command-line arguments:
```bash
# Dry run for a field
casa -c run_imaging.py --field G008.672-00.682 --dry-run

# Run imaging for a field
casa -c run_imaging.py --field G008.672-00.682

# Run imaging for specific spectral window only
casa -c run_imaging.py --field G008.672-00.682 --spw 0
```

### 5. **IMAGING_PLAN_README.md** (4.7 KB)
Detailed documentation of the imaging plan methodology and results.

### 6. **export_regions.py** (8.8 KB)
   - Export imaging groups as DS9 region files
   - Visualize group coverage on astronomical images
   - Supports both DS9 and CASA CRTF formats

### 7. **QUICKSTART.md** (this file)
Quick reference guide.

### 8. **Region Files**
   - `dihca2_imaging_groups.reg` - All 211 groups (25 KB)
   - `dihca2_imaging_groups_with_sources.reg` - Groups + individual sources (68 KB)
   - `region_files/` - Separate region file for each of 65 fields

## Key Statistics

- **Total imaging commands**: 844 (211 groups × 4 SPWs)
- **Total fields**: 65
- **Total sources**: 533
- **Average sources per group**: 2.5
- **Image sizes**:
  - 512×512 pixels (5.12" FOV): 159 groups (75%)
  - 1024×1024 pixels (10.24" FOV): 52 groups (25%)

## Spectral Windows

Each group is imaged in 4 spectral windows:
- **SPW 0**: H₂CO 3(0,3)-2(0,2) @ 218.541687 GHz
- **SPW 1**: ¹³CO 2-1 @ 220.398684 GHz
- **SPW 2**: SiO 5-4 @ 217.104984 GHz
- **SPW 3**: C¹⁸O 2-1 @ 219.560358 GHz

## Visualizing the Groups

DS9 region files have been created to visualize the imaging groups:

```bash
# View all groups in DS9
ds9 your_image.fits -region dihca2_imaging_groups.reg

# View groups with individual source positions
ds9 your_image.fits -region dihca2_imaging_groups_with_sources.reg

# View a specific field
ds9 your_image.fits -region region_files/imaging_groups_008_672m00_682.reg
```

### Generate Custom Region Files

```bash
# Export all groups (basic)
python export_regions.py

# Export with source positions
python export_regions.py --show-sources

# Export only specific field
python export_regions.py --field G008.672-00.682 --show-sources

# Export separate files per field
python export_regions.py --separate-fields --show-sources --output-dir region_files

# Export in CASA CRTF format
python export_regions.py --format crtf --output groups.crtf

# Export in both formats
python export_regions.py --format both
```

## Recommended Workflow

### 1. Explore the Plan
```bash
# See overall statistics
python query_imaging_plan.py --stats

# Check your field of interest
python query_imaging_plan.py --field YOUR_FIELD --verbose
```

### 2. Test with Dry Run
```bash
# Test what would be run for one field
casa -c run_imaging.py --field G008.672-00.682 --dry-run
```

### 3. Run Imaging
```bash
# Image one field
casa -c run_imaging.py --field G008.672-00.682

# Or run a single group for testing
casa --nologger -c "execfile('run_imaging.py'); run_group('G008.672-00.682_group1')"
```

### 4. Batch Processing
For processing many fields, create a batch script:

```bash
#!/bin/bash
# batch_imaging.sh

FIELDS=(
    "G008.672-00.682"
    "G008.872-00.492"
    "G009.213-00.202"
    # ... add more fields
)

for field in "${FIELDS[@]}"; do
    echo "Processing $field"
    casa --nologger --nogui -c "execfile('run_imaging.py'); run_field('$field')"
done
```

## Output Products

Each tclean run will produce:
- `{field}_group{N}_spw{M}.image` - cleaned image cube
- `{field}_group{N}_spw{M}.residual` - residual cube
- `{field}_group{N}_spw{M}.psf` - point spread function
- `{field}_group{N}_spw{M}.pb` - primary beam
- `{field}_group{N}_spw{M}.model` - model cube
- `{field}_group{N}_spw{M}.sumwt` - sum of weights

## Troubleshooting

### Measurement Set Not Found
If you get errors about missing MS files, check:
1. MS files exist in `/orange/adamginsburg/ALMA_IMF/Uploads/dihca2/original_uvdata/`
2. Field names match (especially for `_eb1` suffix fields)
3. Update `find_ms_file()` function in `imaging_plan.py` if needed

### Memory Issues
If imaging runs out of memory:
1. Process smaller fields first
2. Process one SPW at a time using `--spw` flag
3. Reduce `niter` or adjust other parameters in the JSON

### Customize Parameters
Edit the generated JSON file to adjust parameters like:
- `niter` (number of iterations)
- `threshold` (stopping threshold)
- `robust` (Briggs weighting parameter)
- `imsize` (image size)

## Top Fields to Image

Based on number of sources:

1. **G345.493+01.469**: 54 sources in 6 groups
2. **G351.444+00.659**: 40 sources in 11 groups
3. **G034.401+00.226**: 35 sources in 8 groups
4. **G333.018+00.766**: 26 sources in 9 groups
5. **G034.258+00.166**: 19 sources in 6 groups

## Running on SLURM (HiPerGator)

For production imaging on HiPerGator, use the SLURM job runner:

```bash
# Test with dry run
python imaging_job_runner.py --dry-run --max-jobs 5

# Submit next 10 jobs (default - recommended for incremental processing)
python imaging_job_runner.py

# Submit next 20 jobs
python imaging_job_runner.py --max-jobs 20

# Submit all remaining jobs (unlimited)
python imaging_job_runner.py --max-jobs 0

# Submit test jobs for one field
python imaging_job_runner.py --field-filter G008.672-00.682 --max-jobs 4

# Submit only specific SPW
python imaging_job_runner.py --spw-filter 0 --max-jobs 20

# Redo failed jobs
python imaging_job_runner.py --redo-failed

# Validate products and delete blank/failed images
python imaging_job_runner.py --validate-products

# Check job status
squeue -u $USER | grep imaging_

# Cleanup working directory when done
python imaging_job_runner.py --run-cleanup
```

**Key directories:**
- **Working**: `/red/adamginsburg/dihca/workdir_grouped` (fast I/O)
- **Products**: `/orange/adamginsburg/salt/dihca2imaging/grouped_imaging_products` (long-term)
- **Logs**: `/red/adamginsburg/logs`

See `IMAGING_JOBS_README.md` for detailed SLURM documentation.

## Questions?

See detailed documentation:
- `IMAGING_PLAN_README.md` - Methodology and imaging plan details
- `IMAGING_JOBS_README.md` - SLURM job runner and production imaging
- `region_files/README.md` - DS9 region files

For issues or questions, check:
- The log files in `imaging_logs/` (interactive) or `/red/adamginsburg/logs/` (SLURM)
- The metadata in each command for source lists and FOV info
- Use `query_imaging_plan.py` to explore the plan structure
