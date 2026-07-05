#!/usr/bin/env python3
"""
Memory-efficient source cataloging using chunked processing
"""

import os
import glob
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.stats import mad_std
from astropy.wcs import WCS
import astropy.units as u
from scipy.ndimage import label
import warnings

def extract_source_name(filename):
    """Extract source name from filename"""
    basename = os.path.basename(filename)
    source_name = basename.replace('.quick_cont.pbclean.image.fits', '')
    return source_name

def process_fits_file_chunked(filepath, chunk_size=2000):
    """Process a single FITS file using chunked processing"""
    print(f"Processing {os.path.basename(filepath)}...", flush=True)

    # Load the FITS file header first
    with fits.open(filepath) as hdul:
        header = hdul[0].header
        data_shape = hdul[0].data.shape

        # Handle different data shapes
        if len(data_shape) == 4:
            ny, nx = data_shape[2], data_shape[3]
            data_slice = (0, 0, slice(None), slice(None))
        elif len(data_shape) == 3:
            ny, nx = data_shape[1], data_shape[2]
            data_slice = (0, slice(None), slice(None))
        else:
            ny, nx = data_shape[0], data_shape[1]
            data_slice = (slice(None), slice(None))

        print(f"  Image dimensions: {ny} x {nx}", flush=True)

        # First pass: estimate noise from a smaller central region
        center_y, center_x = ny // 2, nx // 2
        noise_region_size = min(1000, ny // 4, nx // 4)
        y1 = center_y - noise_region_size // 2
        y2 = center_y + noise_region_size // 2
        x1 = center_x - noise_region_size // 2
        x2 = center_x + noise_region_size // 2

        if len(data_shape) == 4:
            noise_data = hdul[0].data[0, 0, y1:y2, x1:x2]
        elif len(data_shape) == 3:
            noise_data = hdul[0].data[0, y1:y2, x1:x2]
        else:
            noise_data = hdul[0].data[y1:y2, x1:x2]

        finite_data = noise_data[np.isfinite(noise_data)]

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

        sources = []
        source_name = extract_source_name(filepath)
        source_count = 0

        # Process image in chunks
        n_chunks_y = (ny + chunk_size - 1) // chunk_size
        n_chunks_x = (nx + chunk_size - 1) // chunk_size

        print(f"  Processing in {n_chunks_y} x {n_chunks_x} chunks of size {chunk_size}", flush=True)

        for chunk_y in range(n_chunks_y):
            for chunk_x in range(n_chunks_x):
                # Calculate chunk boundaries
                y_start = chunk_y * chunk_size
                y_end = min((chunk_y + 1) * chunk_size, ny)
                x_start = chunk_x * chunk_size
                x_end = min((chunk_x + 1) * chunk_size, nx)

                # Load chunk data
                if len(data_shape) == 4:
                    chunk_data = hdul[0].data[0, 0, y_start:y_end, x_start:x_end]
                elif len(data_shape) == 3:
                    chunk_data = hdul[0].data[0, y_start:y_end, x_start:x_end]
                else:
                    chunk_data = hdul[0].data[y_start:y_end, x_start:x_end]

                # Process chunk
                clean_data = np.where(np.isfinite(chunk_data), chunk_data, 0)
                mask = clean_data > detection_threshold

                if np.any(mask):
                    # Label connected regions in this chunk
                    labeled_array, num_features = label(mask)

                    # Process each labeled region
                    for i in range(1, num_features + 1):
                        region_mask = labeled_array == i
                        region_data = clean_data[region_mask]

                        # Get peak value in this region
                        peak_value = np.max(region_data)

                        # Only keep sources with peak > 10-sigma
                        if peak_value >= peak_threshold:
                            # Get peak position in chunk coordinates
                            peak_indices = np.where((clean_data == peak_value) & region_mask)
                            if len(peak_indices[0]) > 0:
                                chunk_peak_y, chunk_peak_x = peak_indices[0][0], peak_indices[1][0]

                                # Convert to global image coordinates
                                global_peak_y = y_start + chunk_peak_y
                                global_peak_x = x_start + chunk_peak_x

                                # Convert pixel coordinates to world coordinates
                                try:
                                    sky_coord = wcs.pixel_to_world(global_peak_x, global_peak_y)
                                    ra = sky_coord.ra.deg
                                    dec = sky_coord.dec.deg
                                except Exception as e:
                                    print(f"    Warning: Could not convert coordinates: {e}", flush=True)
                                    ra = np.nan
                                    dec = np.nan

                                # Calculate signal-to-noise ratio
                                snr = peak_value / noise_level

                                # Count pixels in this source
                                npix = np.sum(region_mask)

                                source_count += 1
                                sources.append({
                                    'source_field': source_name,
                                    'source_number': source_count,
                                    'source_id': f"{source_name}_{source_count:03d}",
                                    'ra_deg': ra,
                                    'dec_deg': dec,
                                    'peak_flux': peak_value,
                                    'snr': snr,
                                    'noise_level': noise_level,
                                    'peak_x_pix': global_peak_x,
                                    'peak_y_pix': global_peak_y,
                                    'npix': npix
                                })

                # Progress indicator
                if (chunk_y * n_chunks_x + chunk_x + 1) % 10 == 0:
                    progress = (chunk_y * n_chunks_x + chunk_x + 1) / (n_chunks_y * n_chunks_x) * 100
                    print(f"    Chunk progress: {progress:.1f}% ({len(sources)} sources so far)", flush=True)

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
        sources = process_fits_file_chunked(fits_file)
        all_sources.extend(sources)
        processed_count += 1

        print(f"  Total sources so far: {len(all_sources)}", flush=True)

        # Save intermediate results every 10 files
        if processed_count % 10 == 0:
            print(f"  Completed {processed_count}/{len(fits_files)} files", flush=True)

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
    table.meta['description'] = 'Source catalog from DIHCA continuum images (chunked processing)'
    table.meta['detection_threshold'] = '3-sigma'
    table.meta['peak_threshold'] = '10-sigma'
    table.meta['noise_estimator'] = 'mad_std'
    table.meta['method'] = 'chunked scipy.ndimage connected component labeling'

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
