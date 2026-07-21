# DIHCA2 Imaging Pipeline

This pipeline provides automated imaging for all spectral windows in all MS files from the DIHCA2 dataset using the ACES parallel tclean approach.

## Overview

The pipeline consists of three main components:

1. **MS Inspection** (`inspect_ms_metadata.py`): Examines MS files to extract metadata
2. **Job Runner** (`dihca2_job_runner.py`): Manages SLURM jobs for parallel imaging
3. **Orchestration** (`run_dihca2_imaging.py`): Coordinates the entire process

## Data Location

The pipeline processes MS files located in:
```
/orange/adamginsburg/ALMA_IMF/Uploads/dihca2/original_uvdata/
```

## Quick Start

### Step 1: Inspect MS Files
This must be run in a CASA environment to access MS metadata:

```bash
casa -c "execfile('run_dihca2_imaging.py'); inspect_ms_files()"
```

This will:
- Examine all `.ms` files in the data directory
- Extract spectral window information, field names, and UV ranges
- Calculate optimal imaging parameters (image size, cell size)
- Generate `dihca2_metadata.json` and `dihca2_tclean_commands.json`

### Step 2: Launch Imaging Jobs
```bash
python run_dihca2_imaging.py --run-imaging
```

Optional flags:
- `--dry-run`: Show what would be done without submitting jobs
- `--redo-completed`: Redo jobs that are already completed
- `--redo-failed`: Redo jobs that failed
- `--spw-filter 25 27 29`: Only process specific spectral windows
- `--field-filter G008 G009`: Only process fields containing these strings
- `--max-jobs 1`: Limit number of jobs submitted (useful for testing)

### Step 3: Monitor Jobs
```bash
python run_dihca2_imaging.py --status
```

### Step 4: Cleanup Large Files
After imaging is complete, remove unnecessary large files:
```bash
python run_dihca2_imaging.py --cleanup
```

This removes `.psf`, `.model`, `.residual`, and `.weight` cubes to save disk space.

## Testing

### Incremental Processing

The pipeline supports incremental processing to avoid reprocessing all MS files:

```bash
# Test metadata generation with just 3 MS files
python run_dihca2_imaging.py --inspect --limit 3

# Force update of all entries (even if already processed with correct image sizes)
python run_dihca2_imaging.py --inspect --force-update --limit 5

# Process remaining MS files incrementally (skips already processed files)
python run_dihca2_imaging.py --inspect
```

### Job Testing

Before running the full pipeline on all MS files, it's recommended to test with a single job:

```bash
# Test with just one job
python run_dihca2_imaging.py --run-imaging --max-jobs 1 --dry-run

# If dry-run looks good, submit the test job
python run_dihca2_imaging.py --run-imaging --max-jobs 1

# Monitor the test job
python run_dihca2_imaging.py --status

# If successful, proceed with more jobs or all jobs
python run_dihca2_imaging.py --run-imaging
```

### Intelligent Job Management

The pipeline includes several smart features to optimize job submission and avoid conflicts:

#### Automatic Output Detection
Before submitting any job, the pipeline checks if all output products already exist:
```bash
# Will automatically skip jobs where all array tasks have completed
python run_dihca2_imaging.py --run-imaging
```

#### Merge Job Conflict Prevention
The pipeline prevents new job submissions when merge jobs are running:
- Automatically detects running merge jobs (jobs ending with `_merge`)
- Blocks new submissions until merge operations complete
- Provides clear status messages about running merge jobs

#### Smart Failed Job Resubmission
```bash
# Resubmit only failed jobs (skips completed and running jobs)
python run_dihca2_imaging.py --run-imaging --redo-failed

# Resubmit completed jobs (useful for testing different parameters)
python run_dihca2_imaging.py --run-imaging --redo-completed
```

**Smart Array Task Resubmission**: When using `--redo-failed`, the pipeline automatically:
- Checks which individual array tasks completed successfully by examining output files
- Only resubmits the specific array tasks that failed or didn't complete
- Skips array tasks that produced valid `.image` files
- Uses optimized SLURM array specifications (e.g., `--array=1,3,5-7,10` for selective tasks)

This saves significant compute time by avoiding re-running successful array tasks and preventing resource conflicts.

## Directory Structure

```
/orange/adamginsburg/salt/dihca2imaging/
├── run_dihca2_imaging.py          # Main orchestration script
├── inspect_ms_metadata.py         # MS inspection script
├── dihca2_job_runner.py           # SLURM job management
├── dihca2_metadata.json           # MS metadata (generated)
├── dihca2_tclean_commands.json    # tclean parameters (generated)
├── final_products/                # Final image products
└── README.md                      # This file

/red/adamginsburg/dihca/workdir/   # High-performance working directory
├── G008.672-00.682_spw0/          # Individual job working directories
├── G008.672-00.682_spw1/
└── ...                            # (temporary files, cleaned up after completion)

/red/adamginsburg/logs/            # High-performance log directory
├── casa_log_split_*.log           # CASA split job logs
├── casa_log_line_*.log            # CASA imaging job logs
├── casa_log_merge_*.log           # CASA merge job logs
└── dihca2_*.log                   # SLURM job logs
```

## Pipeline Details

### MS Inspection

The inspection script (`inspect_ms_metadata.py`) extracts:

- **Spectral Windows**: Number of channels, frequency range, channel width
- **Field Information**: Names, coordinates (RA/Dec)
- **UV Coverage**: Min/max baselines for determining image parameters
- **Imaging Parameters**: Optimal cell size and image size based on UV coverage

### Parallel Imaging

The pipeline uses ACES `parallel_tclean` which:

- Splits imaging into channel chunks for parallel processing
- Uses SLURM job arrays for efficient resource utilization
- Automatically merges channel chunks into final cubes
- Handles job dependencies and error recovery
- Uses `/red/adamginsburg/dihca/workdir/` for high-performance I/O during processing
- Moves final products to `/orange/adamginsburg/salt/dihca2imaging/final_products/` for long-term storage

### Default Parameters

```python
default_parameters = {
    'mem': 128,              # GB memory per job
    'ntasks': 16,            # SLURM tasks per job
    'nchan_per': 128,        # Channels per parallel chunk
    'mem_per_cpu': '8gb',    # Memory per CPU
    'jobtime': '96:00:00'    # Maximum job time
}
```

### Image Size Constraints

To manage computational resources, the pipeline enforces:
- **Maximum image size**: 4096×4096 pixels
- **Adaptive cell size**: If the optimal image size exceeds 4096 pixels, the cell size is increased by up to 3× to keep the image size manageable
- **Power-of-2 sizing**: Image sizes are rounded to the nearest power of 2 for computational efficiency

### tclean Parameters

Default tclean parameters include:
- `gridder`: 'mosaic'
- `deconvolver`: 'hogbom'
- `weighting`: 'briggs'
- `robust`: 0.5
- `niter`: 10000
- `threshold`: '0.1mJy'

## File Outputs

For each MS and spectral window, the pipeline produces:

**Essential Files** (kept):
- `.image`: Final cleaned image
- `.image.pbcor`: Primary beam corrected image
- `.pb`: Primary beam response
- `.image.fits`: FITS version of image
- `.image.pbcor.fits`: FITS version of PB-corrected image

**Large Files** (removed by cleanup):
- `.psf`: Point spread function
- `.model`: Model components
- `.residual`: Residual image
- `.weight`: Weighting function

## Monitoring and Troubleshooting

### Check Job Status
```bash
# Check DIHCA2 jobs specifically
python run_dihca2_imaging.py --status

# Check all jobs
squeue -u $USER

# Check job details
scontrol show job <jobid>
```

### Log Files
- SLURM logs: `/red/adamginsburg/logs/dihca2_<field>_spw<N>_*.log`
- CASA logs: `/red/adamginsburg/logs/casa_log_*_<timestamp>.log`

### Common Issues

1. **Out of Memory**: Reduce `nchan_per` parameter
2. **Timeout**: Increase `jobtime` or reduce `nchan_per`
3. **Failed Jobs**: Check logs and use `--redo-failed`

### Disk Space Management

The pipeline can generate large amounts of data. Monitor disk usage:
```bash
du -sh /red/adamginsburg/dihca/workdir/
du -sh /red/adamginsburg/logs/
```

Run cleanup regularly:
```bash
python run_dihca2_imaging.py --cleanup
```

## Dependencies

- CASA (6.4.3-2-pipeline-2021.3.0.17 or compatible)
- ACES reduction pipeline
- SLURM workload manager
- Python 3.6+
- NumPy, Astropy

## SLURM Configuration

The pipeline uses these SLURM settings:
- Account: `astronomy-dept`
- QOS: `astronomy-dept-b`
- Partition: Determined by QOS
- Time limit: 96 hours
- Memory: 8GB per CPU

## Customization

### Modify Imaging Parameters

Edit the tclean parameters in `inspect_ms_metadata.py`:

```python
base_params = {
    'imsize': [metadata['suggested_imaging_params']['suggested_imsize_pixels']] * 2,
    'cell': f"{metadata['suggested_imaging_params']['suggested_cell_arcsec']:.3f}arcsec",
    'weighting': 'briggs',
    'robust': 0.5,
    'niter': 10000,
    'threshold': '0.1mJy',
    # Add or modify parameters here
}
```

### Adjust Job Parameters

Modify `default_parameters` in `dihca2_job_runner.py`:

```python
self.default_parameters = {
    'mem': 256,           # Increase memory
    'ntasks': 32,         # More parallel tasks
    'nchan_per': 64,      # Smaller chunks
    'mem_per_cpu': '8gb',
    'jobtime': '200:00:00'  # Longer time limit
}
```

## Contact

For questions or issues, contact the ACES team or check the ACES documentation.
