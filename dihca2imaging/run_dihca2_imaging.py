#!/usr/bin/env python3
"""
Main orchestration script for DIHCA2 imaging pipeline.

This script coordinates the entire imaging process:
1. Inspects MS files to extract metadata
2. Generates tclean parameters in ACES format
3. Launches parallel imaging jobs
4. Provides cleanup utilities

Usage:
  # Step 1: Inspect MS files (run in CASA)
  casa -c "execfile('run_dihca2_imaging.py'); inspect_ms_files()"

  # Step 2: Launch imaging jobs
  python run_dihca2_imaging.py --run-imaging

  # Step 3: Clean up large files after completion
  python run_dihca2_imaging.py --cleanup
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path


def setup_environment():
    """Set up the environment and check dependencies."""

    # Check that ACES is available
    aces_path = '/orange/adamginsburg/ACES/reduction_ACES'
    if not os.path.exists(aces_path):
        print(f"Error: ACES not found at {aces_path}")
        return False

    # Add ACES to Python path
    if aces_path not in sys.path:
        sys.path.insert(0, aces_path)

    # Check data directory
    data_dir = '/orange/adamginsburg/ALMA_IMF/Uploads/dihca2/original_uvdata'
    if not os.path.exists(data_dir):
        print(f"Error: Data directory not found: {data_dir}")
        return False

    # Create working directories
    base_dir = '/orange/adamginsburg/salt/dihca2imaging'
    high_perf_workdir = '/red/adamginsburg/dihca/workdir'
    high_perf_logdir = '/red/adamginsburg/logs'

    # Create directories on both filesystems
    os.makedirs(os.path.join(base_dir, 'final_products'), exist_ok=True)
    os.makedirs(high_perf_workdir, exist_ok=True)
    os.makedirs(high_perf_logdir, exist_ok=True)

    return True


def inspect_ms_files(limit=None, limit_incomplete=None, force_update=False):
    """
    Run MS inspection to generate metadata.
    This function should be called from within CASA.
    """
    from casatools import msmetadata

    # Import and run the inspection script as subprocess to avoid argument parser conflicts
    inspect_script = '/orange/adamginsburg/salt/dihca2imaging/inspect_ms_metadata.py'

    if not os.path.exists(inspect_script):
        print(f"Error: Inspection script not found: {inspect_script}")
        return False

    # Run the inspection script with optional arguments
    cmd = ['/red/adamginsburg/miniconda3/envs/python312/bin/python', inspect_script]

    if limit:
        cmd.extend(['--limit', str(limit)])
    if limit_incomplete:
        cmd.extend(['--limit-incomplete', str(limit_incomplete)])
    if force_update:
        cmd.append('--force-update')

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        print("MS inspection completed successfully")
        print(result.stdout)
        return True
    else:
        print(f"MS inspection failed with return code {result.returncode}")
        print(f"Error: {result.stderr}")
        return False


def run_imaging_jobs(args):
    """Launch the imaging jobs using the job runner."""

    job_runner_script = '/orange/adamginsburg/salt/dihca2imaging/dihca2_job_runner.py'

    if not os.path.exists(job_runner_script):
        print(f"Error: Job runner script not found: {job_runner_script}")
        return False

    # Check that metadata files exist
    metadata_file = '/orange/adamginsburg/salt/dihca2imaging/dihca2_metadata.json'
    tclean_file = '/orange/adamginsburg/salt/dihca2imaging/dihca2_tclean_commands.json'

    if not os.path.exists(metadata_file):
        print(f"Error: Metadata file not found: {metadata_file}")
        print("Run MS inspection first: casa -c \"execfile('run_dihca2_imaging.py'); inspect_ms_files()\"")
        return False

    if not os.path.exists(tclean_file):
        print(f"Error: Tclean commands file not found: {tclean_file}")
        print("Run MS inspection first: casa -c \"execfile('run_dihca2_imaging.py'); inspect_ms_files()\"")
        return False

    # Build command
    cmd = ['python', job_runner_script]

    if args.dry_run:
        cmd.append('--dry-run')
    if args.redo_completed:
        cmd.append('--redo-completed')
    if args.redo_failed:
        cmd.append('--redo-failed')
    if args.spw_filter:
        cmd.extend(['--spw-filter'] + [str(spw) for spw in args.spw_filter])
    if args.field_filter:
        cmd.extend(['--field-filter'] + args.field_filter)
    if args.max_jobs:
        cmd.extend(['--max-jobs', str(args.max_jobs)])

    print(f"Running: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, check=True)
        print("Imaging jobs launched successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error launching imaging jobs: {e}")
        return False


def run_cleanup():
    """Run cleanup to remove large files."""

    job_runner_script = '/orange/adamginsburg/salt/dihca2imaging/dihca2_job_runner.py'

    if not os.path.exists(job_runner_script):
        print(f"Error: Job runner script not found: {job_runner_script}")
        return False

    # First create the cleanup script
    print("Creating cleanup script...")
    cmd_create = ['python', job_runner_script, '--create-cleanup-script']
    try:
        subprocess.run(cmd_create, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error creating cleanup script: {e}")
        return False

    # Then run the cleanup
    print("Running cleanup...")
    cmd_run = ['python', job_runner_script, '--run-cleanup']
    try:
        subprocess.run(cmd_run, check=True)
        print("Cleanup completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running cleanup: {e}")
        return False


def check_job_status():
    """Check the status of running imaging jobs."""

    try:
        # Get SLURM job status
        result = subprocess.run([
            '/opt/slurm/bin/squeue',
            '--format=%j,%T,%R,%M',
            '--name=dihca2*'
        ], capture_output=True, text=True)

        if result.stdout.strip():
            print("Current DIHCA2 imaging jobs:")
            print("Job Name, State, Reason, Time")
            print(result.stdout)
        else:
            print("No DIHCA2 imaging jobs currently running")

        return True

    except subprocess.CalledProcessError as e:
        print(f"Error checking job status: {e}")
        return False


def print_usage():
    """Print usage information."""

    usage = """
DIHCA2 Imaging Pipeline Usage:

1. INSPECT MS FILES (run in CASA):
   casa -c "execfile('run_dihca2_imaging.py'); inspect_ms_files()"

   This will:
   - Examine all .ms files in /orange/adamginsburg/ALMA_IMF/Uploads/dihca2/original_uvdata
   - Extract metadata (spectral windows, field names, UV range)
   - Generate tclean parameters
   - Save results to JSON files

2. LAUNCH IMAGING JOBS:
   python run_dihca2_imaging.py --run-imaging

   Optional flags:
   --dry-run                  : Show what would be done without submitting jobs
   --redo-completed          : Redo jobs that are already completed
   --redo-failed            : Redo jobs that failed
   --spw-filter 25 27 29    : Only process specific spectral windows
   --field-filter G008 G009 : Only process fields containing these strings
   --max-jobs 1             : Limit number of jobs (useful for testing)

3. CHECK JOB STATUS:
   python run_dihca2_imaging.py --status

4. CLEANUP LARGE FILES:
   python run_dihca2_imaging.py --cleanup

   This removes .psf, .model, .residual, and .weight cubes to save disk space

FILES CREATED:
- dihca2_metadata.json         : MS file metadata
- dihca2_tclean_commands.json  : tclean parameters for each MS/spw
- workdir/                     : Working directory for imaging
- logs/                        : SLURM and CASA log files
- final_products/              : Final image products
"""

    print(usage)


def main():
    parser = argparse.ArgumentParser(
        description='DIHCA2 Imaging Pipeline Orchestration',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Main actions
    parser.add_argument('--inspect', action='store_true',
                        help='Inspect MS files (must be run in CASA)')
    parser.add_argument('--run-imaging', action='store_true',
                        help='Launch parallel imaging jobs')
    parser.add_argument('--cleanup', action='store_true',
                        help='Clean up large files')
    parser.add_argument('--status', action='store_true',
                        help='Check job status')
    parser.add_argument('--usage', action='store_true',
                        help='Show detailed usage information')

    # Job control options
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be done without submitting jobs')
    parser.add_argument('--redo-completed', action='store_true',
                        help='Redo jobs that are already completed')
    parser.add_argument('--redo-failed', action='store_true',
                        help='Redo jobs that failed')

    # Filters
    parser.add_argument('--spw-filter', type=int, nargs='+',
                        help='Only process these spectral windows')
    parser.add_argument('--field-filter', nargs='+',
                        help='Only process fields containing these strings')
    parser.add_argument('--max-jobs', type=int,
                        help='Maximum number of jobs to submit (useful for testing)')
    parser.add_argument('--limit', type=int,
                        help='Limit MS inspection to first N files (for testing)')
    parser.add_argument('--limit-incomplete', type=int,
                        help='Limit MS inspection to next N incomplete files (missing or needing size constraints)')
    parser.add_argument('--force-update', action='store_true',
                        help='Force update of all MS metadata entries')

    args = parser.parse_args()

    # Show usage if requested or no action specified
    if args.usage or not any([args.inspect, args.run_imaging, args.cleanup, args.status]):
        print_usage()
        return

    # Set up environment
    if not setup_environment():
        sys.exit(1)

    # Execute requested actions
    success = True

    if args.inspect:
        success = inspect_ms_files(limit=args.limit, limit_incomplete=args.limit_incomplete, force_update=args.force_update)

    if args.run_imaging:
        success = success and run_imaging_jobs(args)

    if args.cleanup:
        success = success and run_cleanup()

    if args.status:
        success = success and check_job_status()

    if not success:
        sys.exit(1)

    print("Pipeline step completed successfully")


if __name__ == '__main__':
    main()
