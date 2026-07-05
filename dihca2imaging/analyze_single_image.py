#!/usr/bin/env python3
"""
Simple script to analyze a specific image file.
Usage: python analyze_single_image.py <image_file_path>
"""

import os
import sys
import numpy as np
from astropy.table import Table
from spectral_analysis import SpectralAnalyzer

def analyze_specific_file(image_file, source_ra=None, source_dec=None):
    """Analyze a specific image file"""

    if not os.path.exists(image_file):
        print(f"ERROR: Image file not found: {image_file}")
        return

    # Configuration
    catalog_file = "/orange/adamginsburg/salt/dihca2imaging/dihca_source_catalog.fits"
    output_directory = "/orange/adamginsburg/salt/dihca2imaging/debug_single_image"
    os.makedirs(output_directory, exist_ok=True)

    print(f"=== ANALYZING SINGLE IMAGE FILE ===")
    print(f"Image file: {image_file}")
    print(f"File size: {os.path.getsize(image_file)/(1024*1024):.1f} MB")
    print()

    # Initialize analyzer (we need this for the methods)
    analyzer = SpectralAnalyzer(catalog_file, "/dummy", output_directory)

    # If no coordinates provided, use the first source from catalog
    if source_ra is None or source_dec is None:
        first_source = analyzer.catalog[0]
        source_ra = first_source['ra_deg']
        source_dec = first_source['dec_deg']
        print(f"Using coordinates from catalog: RA={source_ra:.6f}, Dec={source_dec:.6f}")
    else:
        print(f"Using provided coordinates: RA={source_ra:.6f}, Dec={source_dec:.6f}")
    print()

    try:
        # Step 1: Load cube
        print("STEP 1: Loading spectral cube...")
        cube = analyzer.load_spectral_cube(image_file)
        print(f"  ✓ Loaded cube: {cube.shape}")
        print(f"  ✓ Units: {cube.unit}")
        print(f"  ✓ Allow huge operations: {cube.allow_huge_operations}")

        # Step 2: Quick data check
        print("\\nSTEP 2: Checking data quality...")
        sample = cube[:50, 100:150, 100:150].filled_data[:].value
        finite_frac = np.sum(np.isfinite(sample)) / sample.size
        print(f"  ✓ Finite data fraction: {finite_frac:.3f}")
        print(f"  ✓ Data range: {np.nanmin(sample):.2e} to {np.nanmax(sample):.2e}")

        # Step 3: Noise estimation
        print("\\nSTEP 3: Estimating noise...")
        from astropy.stats import mad_std
        noise_sample = cube[:100, :50, :50].filled_data[:].value
        finite_noise = noise_sample[np.isfinite(noise_sample)]
        if len(finite_noise) > 100:
            noise_level = mad_std(finite_noise)
            print(f"  ✓ Noise level: {noise_level:.2e}")
        else:
            noise_level = 1e-6
            print(f"  ⚠ Using fallback noise: {noise_level:.2e}")

        # Step 4: Create mask
        print("\\nSTEP 4: Creating source mask...")
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        source_coord = SkyCoord(ra=source_ra*u.deg, dec=source_dec*u.deg)
        mask = analyzer.create_source_mask(cube, source_coord, noise_level)
        mask_pixels = np.sum(mask)
        print(f"  ✓ Mask created: {mask_pixels} pixels selected")

        if mask_pixels == 0:
            print("  ⚠ WARNING: No pixels in mask! This might indicate:")
            print("    - Source coordinates don't match the image")
            print("    - Noise level is too high")
            print("    - Data quality issues")
            print("Because of the cataloging approach we adopted, this outcome should not be reachable."
            raise ValueError("No pixels in mask!")

        # Step 5: Extract spectrum
        print("\\nSTEP 5: Extracting spectrum...")
        spectrum = analyzer.extract_spectrum(cube, mask)
        print(f"  ✓ Extracted spectrum: {len(spectrum)} channels")

        # Check spectrum quality
        if hasattr(spectrum, 'flux'):
            flux_data = spectrum.flux.value
        else:
            flux_data = np.array(spectrum)

        finite_flux = np.sum(np.isfinite(flux_data))
        print(f"  ✓ Finite flux values: {finite_flux}/{len(flux_data)}")

        if finite_flux > 0:
            print(f"  ✓ Flux range: {np.nanmin(flux_data):.2e} to {np.nanmax(flux_data):.2e}")
            print(f"  ✓ Flux mean: {np.nanmean(flux_data):.2e}")
        else:
            print("  ⚠ WARNING: No finite flux values!")
            return

        # Step 6: Find peaks
        print("\\nSTEP 6: Finding spectral peaks...")
        peaks_catalog = analyzer.find_spectral_peaks(spectrum, noise_level)
        print(f"  ✓ Found {len(peaks_catalog)} peaks above 5σ")

        if len(peaks_catalog) > 0:
            print("  Peak summary:")
            for i, peak in enumerate(peaks_catalog[:5]):
                print(f"    Peak {i+1}: channel {peak.get('channel', 'N/A')}, "
                      f"SNR {peak.get('snr', 'N/A'):.1f}")

        # Step 7: Create a simple plot
        print("\\nSTEP 7: Creating diagnostic plot...")
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

        # Plot spectrum
        if hasattr(spectrum, 'spectral_axis'):
            x_axis = spectrum.spectral_axis.value
            ax1.plot(x_axis, flux_data, 'k-', alpha=0.7)
            ax1.set_xlabel('Spectral Axis')
        else:
            ax1.plot(flux_data, 'k-', alpha=0.7)
            ax1.set_xlabel('Channel')

        ax1.set_ylabel('Flux')
        ax1.set_title(f'Extracted Spectrum - {os.path.basename(image_file)}')
        ax1.grid(True, alpha=0.3)

        # Mark peaks
        for peak in peaks_catalog:
            if hasattr(spectrum, 'spectral_axis'):
                peak_x = peak.get('frequency', peak.get('channel', 0))
            else:
                peak_x = peak.get('channel', 0)
            ax1.axvline(peak_x, color='red', alpha=0.7, linestyle='--')

        # Plot mask
        mask_2d = np.sum(mask, axis=0)
        im = ax2.imshow(mask_2d, origin='lower', cmap='viridis')
        ax2.set_title(f'Source Mask ({mask_pixels} pixels)')
        plt.colorbar(im, ax=ax2)

        plt.tight_layout()

        # Save plot
        plot_file = os.path.join(output_directory, f"debug_{os.path.basename(image_file)}.png")
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"  ✓ Saved diagnostic plot: {plot_file}")

        print("\\n=== ANALYSIS SUCCESSFUL ===")
        print(f"Summary:")
        print(f"  - Cube shape: {cube.shape}")
        print(f"  - Noise level: {noise_level:.2e}")
        print(f"  - Mask pixels: {mask_pixels}")
        print(f"  - Spectrum channels: {len(spectrum)}")
        print(f"  - Spectral peaks: {len(peaks_catalog)}")
        print(f"  - Diagnostic plot: {plot_file}")

    except Exception as e:
        print(f"\\n❌ ERROR during analysis: {e}")
        import traceback
        print("\\nFull traceback:")
        traceback.print_exc()

def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_single_image.py <image_file_path> [ra] [dec]")
        print()
        print("Examples:")
        print("  python analyze_single_image.py /path/to/image.fits")
        print("  python analyze_single_image.py /path/to/image.fits 123.456 -12.345")
        print()
        print("If RA/Dec not provided, will use first source from catalog")
        return

    image_file = sys.argv[1]

    # Optional RA/Dec
    source_ra = None
    source_dec = None
    if len(sys.argv) >= 4:
        try:
            source_ra = float(sys.argv[2])
            source_dec = float(sys.argv[3])
        except ValueError:
            print("ERROR: RA and Dec must be numbers")
            return

    analyze_specific_file(image_file, source_ra, source_dec)

if __name__ == "__main__":
    main()
