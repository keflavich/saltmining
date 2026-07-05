#!/usr/bin/env python3
"""
Update the spectral analysis gallery after running new analyses.
This script should be run whenever new diagnostic images are created.
"""

import os
import sys
import subprocess
from pathlib import Path

def update_gallery():
    """Update the gallery HTML file"""

    # Get the directory containing this script
    script_dir = Path(__file__).parent
    results_dir = script_dir / "spectral_analysis_results"

    if not results_dir.exists():
        print(f"Results directory not found: {results_dir}")
        return False

    # Run the gallery generator
    generate_script = script_dir / "generate_gallery.py"

    if not generate_script.exists():
        print(f"Gallery generator not found: {generate_script}")
        return False

    try:
        # Run the generator
        result = subprocess.run([
            sys.executable,
            str(generate_script),
            str(results_dir)
        ], capture_output=True, text=True, check=True)

        print("Gallery updated successfully!")
        print(result.stdout)

        # Print the gallery URL
        gallery_file = results_dir / "gallery.html"
        print(f"\\nView the gallery at:")
        print(f"file://{gallery_file.absolute()}")

        return True

    except subprocess.CalledProcessError as e:
        print(f"Error updating gallery: {e}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        return False

if __name__ == "__main__":
    success = update_gallery()
    if not success:
        sys.exit(1)
