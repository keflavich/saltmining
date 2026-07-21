#!/usr/bin/env python3
"""
Debug script to run spectral analysis on a single image file.
This helps diagnose issues with the spectral analysis pipeline.
"""

import os
import sys
import numpy as np
from astropy.table import Table
from spectral_analysis import SpectralAnalyzer

def debug_single_file():
    """Run analysis on a single file with detailed debugging output"""

    # Configuration
    catalog_file = "/orange/adamginsburg/salt/dihca2imaging/dihca_source_catalog.fits"
    image_directory = "/orange/adamginsburg/salt/dihca2imaging/grouped_imaging_products"
    output_directory = "/orange/adamginsburg/salt/dihca2imaging/debug_analysis"

    # Create debug output directory
    os.makedirs(output_directory, exist_ok=True)

    print("=== SPECTRAL ANALYSIS DEBUG MODE ===")
    print(f"Catalog file: {catalog_file}")
    print(f"Image directory: {image_directory}")
    print(f"Output directory: {output_directory}")
    print()

    # Initialize analyzer
    print("Initializing spectral analyzer...")
    analyzer = SpectralAnalyzer(catalog_file, image_directory, output_directory)
    print(f"Loaded catalog with {len(analyzer.catalog)} sources")
    print()

    # Get the first source for debugging
    first_source = analyzer.catalog[0]
    print("=== SOURCE INFORMATION ===")
    print(f"Source ID: {first_source['source_id']}")
    print(f"Field: {first_source['source_field']}")
    print(f"RA: {first_source['ra_deg']:.6f}")
    print(f"Dec: {first_source['dec_deg']:.6f}")
    print(f"Peak SNR: {first_source['snr']:.2f}")
    print()

    # Find image files for this source
    print("=== FINDING IMAGE FILES ===")
    image_files = analyzer.find_image_files(first_source['source_field'])
    print(f"Found {len(image_files)} image files:")
    for i, img_file in enumerate(image_files):
        print(f"  {i+1}. {os.path.basename(img_file)}")
        print(f"     Full path: {img_file}")
        print(f"     Exists: {os.path.exists(img_file)}")
        if os.path.exists(img_file):
            size_mb = os.path.getsize(img_file) / (1024*1024)
            print(f"     Size: {size_mb:.1f} MB")
    print()

    if not image_files:
        print("ERROR: No image files found!")
        return

    # Select first image file for detailed analysis
    test_image = image_files[0]
    print(f"=== ANALYZING SINGLE FILE ===")
    print(f"Test image: {os.path.basename(test_image)}")
    print()

    try:
        # Load the spectral cube
        print("1. Loading spectral cube...")
        cube = analyzer.load_spectral_cube(test_image)
        print(f"   Cube shape: {cube.shape}")
        print(f"   Cube unit: {cube.unit}")
        print(f"   Has beam: {hasattr(cube, 'beam')}")
        if hasattr(cube, 'beam'):
            print(f"   Beam: {cube.beam}")
        print(f"   Allow huge operations: {cube.allow_huge_operations}")
        print()

        # Sample some data to check for issues
        print("2. Checking cube data...")
        try:
            sample_data = cube[:10, :10, :10].filled_data[:].value
            finite_count = np.sum(np.isfinite(sample_data))
            total_count = sample_data.size
            print(f"   Sample data shape: {sample_data.shape}")
            print(f"   Finite values: {finite_count}/{total_count} ({100*finite_count/total_count:.1f}%)")
            print(f"   Data range: {np.nanmin(sample_data):.2e} to {np.nanmax(sample_data):.2e}")
            print(f"   Data mean: {np.nanmean(sample_data):.2e}")
            print(f"   Data std: {np.nanstd(sample_data):.2e}")
        except Exception as e:
            print(f"   Error sampling data: {e}")
        print()

        # Estimate noise
        print("3. Estimating noise...")
        try:
            # Try to get a small sample from the cube for noise estimation
            sample_data = cube[:100, :50, :50].filled_data[:].value
            finite_sample = sample_data[np.isfinite(sample_data)]
            if len(finite_sample) > 100:
                from astropy.stats import mad_std
                noise_level = mad_std(finite_sample)
                print(f"   Noise level (MAD): {noise_level:.2e}")
            else:
                print(f"   Not enough finite data for noise estimation ({len(finite_sample)} values)")
                noise_level = 1e-6
                print(f"   Using fallback noise: {noise_level:.2e}")
        except Exception as e:
            print(f"   Error estimating noise: {e}")
            noise_level = 1e-6
        print()

        # Create source mask
        print("4. Creating source mask...")
        try:
            from astropy.coordinates import SkyCoord
            import astropy.units as u
            source_coord = SkyCoord(ra=first_source['ra_deg']*u.deg, dec=first_source['dec_deg']*u.deg)
            source_mask = analyzer.create_source_mask(cube, source_coord, noise_level)
            mask_pixels = np.sum(source_mask)
            total_pixels = source_mask.size
            print(f"   Mask shape: {source_mask.shape}")
            print(f"   Mask pixels: {mask_pixels}/{total_pixels} ({100*mask_pixels/total_pixels:.3f}%)")
        except Exception as e:
            print(f"   Error creating mask: {e}")
            source_mask = None
        print()

        # Extract spectrum
        print("5. Extracting spectrum...")
        if source_mask is not None:
            try:
                spectrum = analyzer.extract_spectrum(cube, source_mask)
                print(f"   Spectrum type: {type(spectrum)}")
                if hasattr(spectrum, '__len__'):
                    print(f"   Spectrum length: {len(spectrum)}")
                    if hasattr(spectrum, 'flux'):
                        finite_flux = np.sum(np.isfinite(spectrum.flux))
                        print(f"   Finite flux values: {finite_flux}/{len(spectrum.flux)}")
                        print(f"   Flux range: {np.nanmin(spectrum.flux):.2e} to {np.nanmax(spectrum.flux):.2e}")
                    elif hasattr(spectrum, '__array__'):
                        arr = np.array(spectrum)
                        finite_vals = np.sum(np.isfinite(arr))
                        print(f"   Finite values: {finite_vals}/{len(arr)}")
                        print(f"   Value range: {np.nanmin(arr):.2e} to {np.nanmax(arr):.2e}")
            except Exception as e:
                print(f"   Error extracting spectrum: {e}")
                spectrum = None
        else:
            print("   Skipping spectrum extraction (no mask)")
            spectrum = None
        print()

        # Find spectral peaks
        print("6. Finding spectral peaks...")
        if spectrum is not None:
            try:
                peaks_catalog = analyzer.find_spectral_peaks(spectrum, noise_level)
                print(f"   Found {len(peaks_catalog)} peaks")
                if len(peaks_catalog) > 0:
                    print("   Peak details:")
                    for i, peak in enumerate(peaks_catalog[:3]):  # Show first 3 peaks
                        print(f"     Peak {i+1}: channel={peak.get('channel', 'N/A')}, "
                              f"intensity={peak.get('intensity', 'N/A'):.2e}, "
                              f"SNR={peak.get('snr', 'N/A'):.1f}")
            except Exception as e:
                print(f"   Error finding peaks: {e}")
                peaks_catalog = Table()
        else:
            print("   Skipping peak finding (no spectrum)")
            peaks_catalog = Table()
        print()

        print("=== ANALYSIS COMPLETE ===")
        print("This debug run helps identify where issues occur in the pipeline.")
        print("Check the output above for any errors or unexpected values.")

    except Exception as e:
        print(f"CRITICAL ERROR during analysis: {e}")
        import traceback
        print("Full traceback:")
        traceback.print_exc()

def main():
    """Main function with command line options"""
    if len(sys.argv) > 1:
        if sys.argv[1] == "--help" or sys.argv[1] == "-h":
            print("Usage: python debug_single_file.py [--list-files]")
            print()
            print("Options:")
            print("  --list-files    List all available image files")
            print("  --help, -h      Show this help message")
            return
        elif sys.argv[1] == "--list-files":
            # List all available image files
            image_directory = "/orange/adamginsburg/salt/dihca2imaging/grouped_imaging_products"
            print("Available image files:")
            import glob
            pattern = os.path.join(image_directory, "*.fits")
            files = sorted(glob.glob(pattern))
            for i, f in enumerate(files[:20]):  # Show first 20
                print(f"  {i+1:2d}. {os.path.basename(f)}")
            if len(files) > 20:
                print(f"     ... and {len(files)-20} more files")
            return

    debug_single_file()

if __name__ == "__main__":
    main()
