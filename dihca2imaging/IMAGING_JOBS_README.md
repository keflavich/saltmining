# DIHCA2 Grouped Source Imaging - SLURM Job Runner

## Overview

The `imaging_job_runner.py` script submits CASA imaging jobs to SLURM for the grouped source imaging plan. It uses the same infrastructure as `dihca2_job_runner.py` but is adapted for the new group-based imaging approach.

## Key Features

### High-Performance I/O
- **Working directory**: `/red/adamginsburg/dihca/workdir_grouped` (fast I/O)
- **Log directory**: `/red/adamginsburg/logs` (fast I/O)
- **Final products**: `/orange/adamginsburg/salt/dihca2imaging/grouped_imaging_products` (long-term storage)

Products are automatically moved from `/red` to `/orange` when imaging completes.

### Memory Management
- **Small images** (512×512): 64GB RAM, 8 tasks
- **Large images** (1024×1024): 96GB RAM, 12 tasks
- Memory per CPU: 8GB

### Job Array Structure
Each imaging job consists of three stages:
1. **Split job**: Splits MS by spectral window (if needed)
2. **Array job**: Parallel imaging across channel chunks (240 tasks for 3840 channels with 16 channels per task)
3. **Merge job**: Combines channel chunks and exports to `/orange`

### Intelligent Job Management
- Checks SLURM queue before submission
- Skips already-completed jobs
- Selective resubmission of failed array tasks
- Prevents duplicate submissions (checks for jobs with any suffix: _split, _arr, _merge)

## Usage

### Basic Usage

```bash
# Test with dry run
python imaging_job_runner.py --dry-run --max-jobs 5

# Submit next 10 jobs (default behavior)
python imaging_job_runner.py

# Submit next 20 jobs
python imaging_job_runner.py --max-jobs 20

# Submit all remaining jobs (unlimited)
python imaging_job_runner.py --max-jobs 0
```

**Default behavior:** By default, only the next 10 jobs that need to be run will be submitted. This prevents overwhelming the queue and allows for incremental processing. Jobs are automatically skipped if:
- Working directory already contains output files
- Job is already running or pending in SLURM (checks for job name with any suffix: `_split`, `_arr`, `_merge`)
- All output products already exist in savedir

**Note:** SLURM jobs are created with suffixes (e.g., `imaging_G009.879-00.751_group2_spw3_arr` for array jobs). The system automatically detects and skips jobs with these suffixes.

### Filter by Field

```bash
# Image only a specific field
python imaging_job_runner.py --field-filter G008.672-00.682

# Image multiple fields
python imaging_job_runner.py --field-filter G008 G009 G010
```

### Filter by Spectral Window

```bash
# Image only SPW 0 (H2CO)
python imaging_job_runner.py --spw-filter 0

# Image only SPW 0 and 1
python imaging_job_runner.py --spw-filter 0 1

# Combine with field filter and limit
python imaging_job_runner.py --field-filter G008 --spw-filter 0 --max-jobs 10
```

### Handling Failed Jobs

```bash
# Redo only failed jobs (automatically skips completed array tasks)
python imaging_job_runner.py --redo-failed

# Redo everything including completed jobs
python imaging_job_runner.py --redo-completed
```

### Cleanup

```bash
# Create cleanup script
python imaging_job_runner.py --create-cleanup-script

# Run cleanup to remove large intermediate files from /red
python imaging_job_runner.py --run-cleanup
```

### Validate Products

The job runner automatically detects and handles blank/failed images:

```bash
# Validate all existing products (checks for blank images and deletes them)
python imaging_job_runner.py --validate-products

# Validate specific field only
python imaging_job_runner.py --validate-products --field-filter G008.672

# Validate specific SPW only
python imaging_job_runner.py --validate-products --spw-filter 0
```

**Blank image detection:**
- Automatically runs when checking completed array tasks
- Uses `spectral-cube` for memory-efficient validation
- Checks file size before validation (max 8 GB in memory)
- For large files (>8 GB), samples first/middle/last channels
- Very large files (>16 GB) skip detailed validation
- Checks for: all zeros, all NaNs, no variation, suspiciously low RMS
- Deletes blank images and marks tasks for reprocessing
- Also validates final merged products
- Handles MemoryError gracefully (skips validation rather than crashing)

## Directory Structure

### Working Directory (/red - temporary, high-performance)
```
/red/adamginsburg/dihca/workdir_grouped/
├── G008.672-00.682_group1_spw0/
│   ├── *.image (channel chunks: 0000.016.image, 0016.016.image, etc.)
│   ├── *.psf
│   ├── *.residual
│   ├── *.model
│   └── *.weight
├── G008.672-00.682_group1_spw1/
└── ...
```

### Final Products (/orange - long-term storage)
```
/orange/adamginsburg/salt/dihca2imaging/grouped_imaging_products/
├── G008.672-00.682_group1_spw0.image (merged cube)
├── G008.672-00.682_group1_spw0.psf
├── G008.672-00.682_group1_spw1.image
└── ...
```

### Logs (/red - temporary)
```
/red/adamginsburg/logs/
├── imaging_G008.672-00.682_group1_spw0_<jobid>_<arrayid>.log
└── ...
```

## Monitoring Jobs

### Check SLURM Queue

```bash
# All jobs
squeue -u $USER

# Only imaging jobs
squeue -u $USER | grep imaging_

# Only merge jobs (critical - blocks new submissions)
squeue -u $USER | grep merge
```

### Check Job Status

```bash
# Recent jobs
sacct --format=JobID,JobName%60,State,Elapsed,CPUTime

# Specific job
sacct -j <jobid> --format=JobID,JobName%60,State,Elapsed,ExitCode
```

### Check Output

```bash
# List completed products
ls -lh /orange/adamginsburg/salt/dihca2imaging/grouped_imaging_products/*.image

# Check log files
tail /red/adamginsburg/logs/imaging_G008.672-00.682_group1_spw0_*.log
```

## Troubleshooting

### Blank/Failed Images

**Symptoms:** Final image products are all zeros, NaNs, or show no variation

**Causes:**
- Imaging failed but SLURM job reported success
- Insufficient data in the measurement set
- Phase center issues
- Corrupt MS file

**Solution:**
```bash
# Validate all products and delete blank ones
python imaging_job_runner.py --validate-products

# After validation, resubmit failed groups
python imaging_job_runner.py --redo-failed
```

The system automatically:
- Detects blank images when checking array tasks
- Deletes blank images and marks for reprocessing
- Validates final merged products
- Reports detailed reasons (all zeros, NaNs, no variation, etc.)

### Job Fails with Memory Error
- Large images may need more memory
- Edit `imaging_job_runner.py` to increase memory limits
- Or reduce `nchan_per` to process fewer channels per task

### Job Stuck in PENDING
- Check queue limits: `squeue -u $USER | wc -l`
- Check QOS limits: `sacctmgr show qos astronomy-dept-b`
- May need to wait for other jobs to complete

### Partial Array Completion
- Run with `--redo-failed` to resubmit only missing tasks
- The script automatically detects which tasks completed
- Blank images are automatically detected and deleted
- Only missing channel chunks will be resubmitted

### Disk Space Issues on /red
- Check usage: `du -sh /red/adamginsburg/dihca/workdir_grouped`
- Run cleanup: `python imaging_job_runner.py --run-cleanup`
- Manually remove old working directories if needed

## Performance Estimates

### Per Group Imaging Time
- **Small groups** (1-2 sources, 512×512): ~24 hours
- **Large groups** (5+ sources, 1024×1024): ~36-48 hours

### Total Project
- 211 groups × 4 SPWs = 844 imaging jobs
- Estimated wall time: 2-4 weeks (depending on queue availability)
- Can run many jobs in parallel to speed up

### Resource Usage
- Working directory space: ~500GB - 1TB (temporary)
- Final products: ~200GB - 500GB (permanent)
- Peak memory per job: 64-96GB

## Advanced Options

### Custom Directories

```bash
python imaging_job_runner.py \
    --workdir /red/custom/workdir \
    --logdir /red/custom/logs \
    --savedir /orange/custom/output
```

### Custom SLURM Settings

```bash
python imaging_job_runner.py \
    --account my-account \
    --qos my-qos \
    --casa-version casa-6.5.0
```

### Custom Imaging Plan

```bash
# Use a modified imaging plan
python imaging_job_runner.py --imaging-plan my_custom_plan.json
```

## Batch Processing Examples

### Submit in Stages

```bash
# Stage 1: Small fields only (faster)
python imaging_job_runner.py --field-filter G008 G009 G010

# Stage 2: Medium fields
python imaging_job_runner.py --field-filter G034 G035

# Stage 3: Large fields
python imaging_job_runner.py --field-filter G345 G351
```

### Process One SPW at a Time

```bash
# Week 1: SPW 0 (H2CO)
python imaging_job_runner.py --spw-filter 0

# Week 2: SPW 1 (13CO)
python imaging_job_runner.py --spw-filter 1

# Week 3: SPW 2 (SiO)
python imaging_job_runner.py --spw-filter 2

# Week 4: SPW 3 (C18O)
python imaging_job_runner.py --spw-filter 3
```

### Test-Production Workflow

```bash
# 1. Test with one group
python imaging_job_runner.py --dry-run --field-filter G008.672-00.682 --spw-filter 0 --max-jobs 1

# 2. Submit test job
python imaging_job_runner.py --field-filter G008.672-00.682 --spw-filter 0 --max-jobs 1

# 3. Monitor until complete (check logs, output products)

# 4. If successful, submit more jobs
python imaging_job_runner.py --spw-filter 0 --max-jobs 50

# 5. Scale up to full production
python imaging_job_runner.py
```

## Comparison with Original Job Runner

| Feature | `dihca2_job_runner.py` | `imaging_job_runner.py` |
|---------|------------------------|-------------------------|
| Input | MS metadata + tclean commands | Grouped imaging plan |
| Fields | Full field imaging | Source groups only |
| Image sizes | 4096×4096 pixels | 512-1536 pixels |
| Memory | 128GB | 64-96GB |
| Job time | 96 hours | 48 hours |
| Working dir | `/red/adamginsburg/dihca/workdir` | `/red/adamginsburg/dihca/workdir_grouped` |
| Output dir | `final_products/` | `grouped_imaging_products/` |

## Related Files

- **imaging_plan.py**: Generates the imaging plan
- **dihca2_imaging_plan.json**: Input imaging plan (844 commands)
- **query_imaging_plan.py**: Query/analyze the imaging plan
- **export_regions.py**: Visualize groups as DS9 regions
- **run_imaging.py**: Run imaging interactively in CASA (not for SLURM)

## Support

For issues or questions:
1. Check log files in `/red/adamginsburg/logs/`
2. Verify input files exist and are valid
3. Check SLURM queue status: `squeue -u $USER`
4. Review this documentation
5. Consult `dihca2_job_runner.py` for similar patterns
