#!/usr/bin/env python3
"""
Simple source cataloging using scipy peak detection instead of astrodendro
"""

import os
import glob
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.stats import mad_std
from astropy.wcs import WCS
import astropy.units as u
from scipy import ndimage
from scipy.ndimage import label
import warnings

def extract_source_name(filename):
    """Extract source name from filename"""
    basename = os.path.basename(filename)
    source_name = basename.replace('.quick_cont.pbclean.image.fits', '')
    return source_name

def process_fits_file_simple(filepath):
    """Process a single FITS file using simple peak detection"""
    print(f"Processing {os.path.basename(filepath)}...", flush=True)

    # Load the FITS file
    with fits.open(filepath) as hdul:
        data = hdul[0].data
        header = hdul[0].header

        # Handle different data shapes
        if data.ndim == 4:
            data = data[0, 0, :, :]
        elif data.ndim == 3:
            data = data[0, :, :]

        # Remove NaN values for noise calculation
        finite_data = data[np.isfinite(data)]

        if len(finite_data) == 0:
            print(f"Warning: No finite data in {filepath}", flush=True)
            return []

        # Estimate noise using mad_std
        noise_level = mad_std(finite_data)
        print(f"  Noise level (mad_std): {noise_level:.6e}", flush=True)

        # Set thresholds
        detection_threshold = 3.0 * noise_level
        peak_threshold = 10.0 * noise_level

        print(f"  Detection threshold (3σ): {detection_threshold:.6e}", flush=True)
        print(f"  Peak threshold (10σ): {peak_threshold:.6e}", flush=True)

        # Create WCS object for coordinate conversion
        wcs = WCS(header).celestial

        # Create a mask for pixels above detection threshold
        clean_data = np.where(np.isfinite(data), data, 0)
        mask = clean_data > detection_threshold

        # Label connected regions
        labeled_array, num_features = label(mask)
        print(f"  Found {num_features} potential source regions", flush=True)

        sources = []
        source_name = extract_source_name(filepath)

        # Process each labeled region
        for i in range(1, num_features + 1):
            region_mask = labeled_array == i
            region_data = clean_data[region_mask]

            # Get peak value in this region
            peak_value = np.max(region_data)

            # Only keep sources with peak > 10-sigma
            if peak_value >= peak_threshold:
                # Get peak position
                peak_indices = np.where((clean_data == peak_value) & region_mask)
                if len(peak_indices[0]) > 0:
                    peak_y, peak_x = peak_indices[0][0], peak_indices[1][0]

                    # Convert pixel coordinates to world coordinates
                    try:
                        sky_coord = wcs.pixel_to_world(peak_x, peak_y)
                        ra = sky_coord.ra.deg
                        dec = sky_coord.dec.deg
                    except Exception as e:
                        print(f"    Warning: Could not convert coordinates for source {i}: {e}", flush=True)
                        ra = np.nan
                        dec = np.nan

                    # Calculate signal-to-noise ratio
                    snr = peak_value / noise_level

                    # Count pixels in this source
                    npix = np.sum(region_mask)

                    sources.append({
                        'source_field': source_name,
                        'source_number': i,
                        'source_id': f"{source_name}_{i:03d}",
                        'ra_deg': ra,
                        'dec_deg': dec,
                        'peak_flux': peak_value,
                        'snr': snr,
                        'noise_level': noise_level,
                        'peak_x_pix': peak_x,
                        'peak_y_pix': peak_y,
                        'npix': npix
                    })

                    print(f"    Source {i}: SNR={snr:.1f}, Peak={peak_value:.6e}, RA={ra:.6f}, Dec={dec:.6f}", flush=True)

        print(f"  Kept {len(sources)} sources above 10σ threshold", flush=True)

        return sources

def main():
    """Process all FITS files"""

    # Find all FITS files
    fits_pattern = "/orange/adamginsburg/salt/dihca/*quick_cont*fits"
    fits_files = glob.glob(fits_pattern)

    print(f"Found {len(fits_files)} FITS files to process", flush=True)

    all_sources = []
    processed_count = 0

    # Process each FITS file
    for i, fits_file in enumerate(sorted(fits_files)):
        print(f"\nProcessing {i+1}/{len(fits_files)}: {os.path.basename(fits_file)}", flush=True)
        sources = process_fits_file_simple(fits_file)
        all_sources.extend(sources)
        processed_count += 1

        # Progress update every 10 files
        if processed_count % 10 == 0:
            print(f"  Processed {processed_count} files so far, found {len(all_sources)} total sources", flush=True)

    if not all_sources:
        print("No sources found!")
        return

    print(f"\nTotal sources found: {len(all_sources)}", flush=True)

    # Create astropy table
    table = Table()

    # Add columns
    table['source_field'] = [s['source_field'] for s in all_sources]
    table['source_number'] = [s['source_number'] for s in all_sources]
    table['source_id'] = [s['source_id'] for s in all_sources]
    table['ra_deg'] = [s['ra_deg'] for s in all_sources] * u.deg
    table['dec_deg'] = [s['dec_deg'] for s in all_sources] * u.deg
    table['peak_flux'] = [s['peak_flux'] for s in all_sources]
    table['snr'] = [s['snr'] for s in all_sources]
    table['noise_level'] = [s['noise_level'] for s in all_sources]
    table['peak_x_pix'] = [s['peak_x_pix'] for s in all_sources]
    table['peak_y_pix'] = [s['peak_y_pix'] for s in all_sources]
    table['npix'] = [s['npix'] for s in all_sources]

    # Add metadata
    table.meta['description'] = 'Source catalog from DIHCA continuum images (simple peak detection)'
    table.meta['detection_threshold'] = '3-sigma'
    table.meta['peak_threshold'] = '10-sigma'
    table.meta['noise_estimator'] = 'mad_std'
    table.meta['method'] = 'scipy.ndimage connected component labeling'

    # Save as FITS table
    output_filename = "/orange/adamginsburg/salt/dihca2imaging/dihca_source_catalog.fits"
    table.write(output_filename, format='fits', overwrite=True)

    print(f"\nCatalog saved to: {output_filename}", flush=True)
    print(f"Catalog contains {len(table)} sources", flush=True)

    # Print summary statistics
    print(f"\nSummary statistics:", flush=True)
    print(f"  SNR range: {np.min(table['snr']):.1f} - {np.max(table['snr']):.1f}")
    print(f"  Peak flux range: {np.min(table['peak_flux']):.6e} - {np.max(table['peak_flux']):.6e}", flush=True)
    print(f"  Number of fields with sources: {len(np.unique(table['source_field']))}", flush=True)
    print(f"  Number of sources: {len(table)}", flush=True)

    # Calculate fields without sources
    all_fields = set()
    for fits_file in sorted(fits_files):
        field_name = extract_source_name(fits_file)
        all_fields.add(field_name)

    fields_with_sources = set(table['source_field'])
    fields_without_sources = all_fields - fields_with_sources

    print(f"  Number of fields without sources: {len(fields_without_sources)}", flush=True)
    if fields_without_sources:
        print(f"  Fields without sources: {', '.join(sorted(fields_without_sources))}", flush=True)

    return table

if __name__ == "__main__":
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    catalog = main()
