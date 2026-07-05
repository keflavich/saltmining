#!/usr/bin/env python
"""
Script to execute CASA tclean imaging commands from the imaging plan.

This script is designed to be run within CASA:
    casa --nologger --nogui -c run_imaging.py

Or with specific field/group filters:
    casa --nologger --nogui -c "execfile('run_imaging.py'); run_field('G008.672-00.682')"

Usage patterns:
    1. Run all imaging commands (careful, this is 844 commands!)
       casa -c run_imaging.py --all

    2. Run imaging for a specific field
       casa -c run_imaging.py --field G008.672-00.682

    3. Run a specific group and spectral window
       casa -c run_imaging.py --group G008.672-00.682_group1 --spw 0

    4. Dry run to see what would be executed
       casa -c run_imaging.py --field G008.672-00.682 --dry-run
"""

import json
import os
import sys
import argparse
from datetime import datetime


def load_imaging_plan(filename='dihca2_imaging_plan.json'):
    """Load the imaging plan JSON file."""
    with open(filename, 'r') as f:
        return json.load(f)


def filter_commands(imaging_plan, field=None, group=None, spw=None):
    """
    Filter imaging commands based on criteria.

    Parameters:
    -----------
    imaging_plan : dict
        Full imaging plan
    field : str, optional
        Field name to filter by
    group : str, optional
        Group name to filter by
    spw : int, optional
        Spectral window to filter by

    Returns:
    --------
    dict : Filtered imaging commands
    """
    filtered = {}

    for key, cmd in imaging_plan.items():
        # Check field
        if field and field.lower() not in key.lower():
            continue

        # Check group
        if group and group.lower() not in key.lower():
            continue

        # Check spw
        if spw is not None and cmd['spw'] != spw:
            continue

        filtered[key] = cmd

    return filtered


def run_tclean(imagename, params, dry_run=False, logfile=None):
    """
    Run CASA tclean with the given parameters.

    Parameters:
    -----------
    imagename : str
        Name for this imaging command (used for logging)
    params : dict
        tclean parameters
    dry_run : bool
        If True, just print what would be executed
    logfile : file object, optional
        Log file to write progress

    Returns:
    --------
    bool : True if successful, False otherwise
    """
    # Remove metadata before passing to tclean
    tclean_params = {k: v for k, v in params.items() if k != 'metadata'}

    """
    If outfile exists, report success / completed
    """

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    msg = f"[{timestamp}] Imaging {imagename}"

    print("\n" + "="*60)
    print(msg)
    print("="*60)

    if logfile:
        logfile.write(msg + "\n")
        logfile.flush()

    if dry_run:
        print("DRY RUN - Would execute tclean with parameters:")
        for key, val in tclean_params.items():
            print(f"  {key}: {val}")
        return True

    try:
        # Check if MS file exists
        vis_file = tclean_params['vis'][0]
        if not os.path.exists(vis_file):
            error_msg = f"ERROR: Measurement set not found: {vis_file}"
            print(error_msg)
            if logfile:
                logfile.write(error_msg + "\n")
                logfile.flush()
            return False

        # Execute tclean
        print(f"Starting tclean for {imagename}...")
        print(f"  Image size: {tclean_params['imsize']}")
        print(f"  Phase center: {tclean_params['phasecenter']}")
        print(f"  SPW: {tclean_params['spw']}")
        print(f"  Number of sources: {params.get('metadata', {}).get('n_sources', 'unknown')}")

        raise NotImplementedError("tclean is not supposed to be run directly....")
        ## Import tclean from CASA
        #from casatasks import tclean

        # # Run tclean
        # tclean(**tclean_params)

        success_msg = f"Successfully completed {imagename}"
        print(success_msg)
        if logfile:
            logfile.write(success_msg + "\n")
            logfile.flush()

        return True

    except Exception as e:
        error_msg = f"ERROR imaging {imagename}: {str(e)}"
        print(error_msg)
        if logfile:
            logfile.write(error_msg + "\n")
            logfile.flush()
        return False


def run_imaging(imaging_plan_file='dihca2_imaging_plan.json',
                field=None, group=None, spw=None,
                dry_run=False, log_dir='imaging_logs'):
    """
    Main function to run imaging from the plan.

    Parameters:
    -----------
    imaging_plan_file : str
        Path to the imaging plan JSON file
    field : str, optional
        Filter by field name
    group : str, optional
        Filter by group name
    spw : int, optional
        Filter by spectral window
    dry_run : bool
        If True, don't actually run tclean
    log_dir : str
        Directory to store log files
    """
    # Load imaging plan
    print(f"Loading imaging plan from {imaging_plan_file}...")
    imaging_plan = load_imaging_plan(imaging_plan_file)
    print(f"Loaded {len(imaging_plan)} total imaging commands")

    # Filter commands
    commands = filter_commands(imaging_plan, field=field, group=group, spw=spw)
    print(f"Filtered to {len(commands)} commands to execute")

    if len(commands) == 0:
        print("No commands match the filter criteria!")
        return

    # Create log directory
    if not dry_run:
        os.makedirs(log_dir, exist_ok=True)
        log_file_path = os.path.join(log_dir,
                                      f"imaging_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
        logfile = open(log_file_path, 'w')
        print(f"Writing log to {log_file_path}")
    else:
        logfile = None
        print("DRY RUN MODE - No actual imaging will be performed")

    # Execute commands
    n_success = 0
    n_failed = 0

    try:
        for i, (key, params) in enumerate(commands.items(), start=1):
            print(f"\n{'='*60}")
            print(f"Command {i}/{len(commands)}: {key}")
            print(f"{'='*60}")

            success = run_tclean(key, params, dry_run=dry_run, logfile=logfile)

            if success:
                n_success += 1
            else:
                n_failed += 1

    finally:
        if logfile:
            logfile.close()

    # Summary
    print("\n" + "="*60)
    print("IMAGING SUMMARY")
    print("="*60)
    print(f"Total commands: {len(commands)}")
    print(f"Successful: {n_success}")
    print(f"Failed: {n_failed}")

    if not dry_run and logfile:
        print(f"Log written to: {log_file_path}")


def run_field(field_name, dry_run=False):
    """Convenience function to run all imaging for a specific field."""
    run_imaging(field=field_name, dry_run=dry_run)


def run_group(group_name, dry_run=False):
    """Convenience function to run all imaging for a specific group."""
    run_imaging(group=group_name, dry_run=dry_run)


if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Run CASA imaging from the imaging plan',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--all', action='store_true',
                        help='Run all imaging commands (844 total!)')
    parser.add_argument('--field', type=str,
                        help='Run imaging for specific field')
    parser.add_argument('--group', type=str,
                        help='Run imaging for specific group')
    parser.add_argument('--spw', type=int,
                        help='Run imaging for specific spectral window')
    parser.add_argument('--dry-run', action='store_true',
                        help='Dry run - show what would be executed without running')
    parser.add_argument('--input', type=str, default='dihca2_imaging_plan.json',
                        help='Input imaging plan JSON file')
    parser.add_argument('--log-dir', type=str, default='imaging_logs',
                        help='Directory for log files')

    args = parser.parse_args()

    # Determine what to run
    if not (args.all or args.field or args.group):
        parser.print_help()
        print("\nERROR: Must specify --all, --field, or --group")
        sys.exit(1)

    if args.all and not args.dry_run:
        # Safety check for running all commands
        response = input("WARNING: This will run 844 imaging commands. Continue? (yes/no): ")
        if response.lower() != 'yes':
            print("Aborted.")
            sys.exit(0)

    # Run imaging
    run_imaging(
        imaging_plan_file=args.input,
        field=args.field,
        group=args.group,
        spw=args.spw,
        dry_run=args.dry_run,
        log_dir=args.log_dir
    )
