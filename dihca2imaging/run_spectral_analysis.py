#!/usr/bin/env python3
"""
Run the full spectral analysis on all DIHCA sources
"""

from spectral_analysis import SpectralAnalyzer
import sys

def main():
    """Run the complete analysis"""

    # Configuration
    catalog_file = "/orange/adamginsburg/salt/dihca2imaging/dihca_source_catalog.fits"
    image_directory = "/orange/adamginsburg/salt/dihca2imaging/grouped_imaging_products"
    output_directory = "/orange/adamginsburg/salt/dihca2imaging/spectral_analysis_results"

    # Initialize analyzer
    print("Initializing spectral analyzer...")
    analyzer = SpectralAnalyzer(catalog_file, image_directory, output_directory)

    # Check command line arguments
    if len(sys.argv) > 1:
        try:
            max_sources = int(sys.argv[1])
            print(f"Running analysis on first {max_sources} sources...")
            analyzer.analyze_all_sources(max_sources=max_sources)
        except ValueError:
            print("Invalid number of sources specified. Using default (all sources).")
            analyzer.analyze_all_sources()
    else:
        print("Running analysis on all sources...")
        analyzer.analyze_all_sources()

if __name__ == "__main__":
    main()
