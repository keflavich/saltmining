#!/usr/bin/env python3
"""
Run the spectral analysis on a larger batch of sources for production
"""

from spectral_analysis import SpectralAnalyzer
import sys
import time

def main():
    """Run the analysis on a batch of sources"""

    # Configuration
    catalog_file = "/orange/adamginsburg/salt/dihca2imaging/dihca_source_catalog.fits"
    image_directory = "/orange/adamginsburg/salt/dihca2imaging/grouped_imaging_products"
    output_directory = "/orange/adamginsburg/salt/dihca2imaging/spectral_analysis_results"

    print(f"Running spectral analysis on all sources...")
    start_time = time.time()

    # Initialize analyzer
    analyzer = SpectralAnalyzer(catalog_file, image_directory, output_directory)

    # Run the analysis
    analyzer.analyze_all_sources()

    end_time = time.time()
    elapsed = end_time - start_time

    print(f"\nAnalysis completed in {elapsed:.1f} seconds")
    print(f"Average time per source: {elapsed/max_sources:.1f} seconds")
    print(f"\nResults available in: {output_directory}")

    # Print some summary statistics
    import os
    result_dirs = [d for d in os.listdir(output_directory)
                   if os.path.isdir(os.path.join(output_directory, d)) and d.startswith('G')]

    print(f"\nProcessed {len(result_dirs)} source directories:")
    for source_dir in result_dirs[:5]:  # Show first 5
        source_path = os.path.join(output_directory, source_dir)
        files = os.listdir(source_path)
        spectrum_files = [f for f in files if 'spectrum.fits' in f]
        diagnostic_files = [f for f in files if 'diagnostics.png' in f]
        print(f"  {source_dir}: {len(spectrum_files)} spectra, {len(diagnostic_files)} diagnostics")

    if len(result_dirs) > 5:
        print(f"  ... and {len(result_dirs) - 5} more")

if __name__ == "__main__":
    main()
