# DIHCA2 Imaging Pipeline

"""
DIHCA2 imaging pipeline using ACES parallel tclean approach.

This module provides automated imaging for all spectral windows in all MS files
from the DIHCA2 dataset located at:
/orange/adamginsburg/ALMA_IMF/Uploads/dihca2/original_uvdata

PIPELINE COMPONENTS:
1. inspect_ms_metadata.py - Inspects MS files to extract metadata
2. dihca2_job_runner.py - Manages SLURM jobs for parallel imaging
3. run_dihca2_imaging.py - Main orchestration script
4. cleanup utilities - Removes large intermediate files

USAGE:
See README.md for detailed usage instructions.

Quick start:
1. casa -c "execfile('run_dihca2_imaging.py'); inspect_ms_files()"
2. python run_dihca2_imaging.py --run-imaging
3. python run_dihca2_imaging.py --cleanup

The pipeline uses ACES parallel_tclean for efficient imaging of large spectral cubes.
"""

import os
import sys

import aces
# Do not add ACES to Python path
# ACES_PATH = '/orange/adamginsburg/ACES/reduction_ACES'
# if ACES_PATH not in sys.path:
#     sys.path.insert(0, ACES_PATH)

# Pipeline configuration
DATA_DIR = '/orange/adamginsburg/ALMA_IMF/Uploads/dihca2/original_uvdata'
PIPELINE_DIR = '/orange/adamginsburg/salt/dihca2imaging'

# Check if pipeline is set up
def check_pipeline_setup():
    """Check if the pipeline is properly set up."""
    required_files = [
        'run_dihca2_imaging.py',
        'inspect_ms_metadata.py',
        'dihca2_job_runner.py',
        'README.md'
    ]

    missing_files = []
    for filename in required_files:
        filepath = os.path.join(PIPELINE_DIR, filename)
        if not os.path.exists(filepath):
            missing_files.append(filename)

    if missing_files:
        print(f"Pipeline setup incomplete. Missing files: {missing_files}")
        return False

    print("Pipeline setup complete")
    return True

if __name__ == '__main__':
    check_pipeline_setup()