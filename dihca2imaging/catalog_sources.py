"""
Catalog all sources in /orange/adamginsburg/salt/dihca/*quick_cont*fits

Use astropy's mad_std to estimate noise, then use astrodendro  with a threshold of 3-sigma to find sources.

Reject all sources that have a peak less than 10-sigma.


Create a catalog of these sources.  They should be labeled by the source name (e.g., G351.444+00.659) and source number in that field.

Save the resulting catalog as an astropy FITS table.
"""

import os
import glob
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
from astropy.stats import mad_std
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astrodendro import Dendrogram
import warnings

def extract_source_name(filename):
    """Extract source name from filename (e.g., G351.444+00.659 from G351.444+00.659.quick_cont.pbclean.image.fits)"""
    basename = os.path.basename(filename)
    # Remove the .quick_cont.pbclean.image.fits suffix
    source_name = basename.replace('.quick_cont.pbclean.image.fits', '')
    return source_name

def process_fits_file(filepath):
    """Process a single FITS file to find sources using astrodendro"""
    print(f"Processing {os.path.basename(filepath)}...")

    # Load the FITS file
    with fits.open(filepath) as hdul:
        data = hdul[0].data
        header = hdul[0].header

        # Handle different data shapes (sometimes 4D, sometimes 2D)
        if data.ndim == 4:
            # Take the first frequency and Stokes parameter
            data = data[0, 0, :, :]
        elif data.ndim == 3:
            # Take the first slice
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
        detection_threshold = 3.0 * noise_level  # 3-sigma for detection
        peak_threshold = 10.0 * noise_level      # 10-sigma for keeping sources

        print(f"  Detection threshold (3σ): {detection_threshold:.6e}", flush=True)
        print(f"  Peak threshold (10σ): {peak_threshold:.6e}", flush=True)

        # Create WCS object for coordinate conversion
        wcs = WCS(header).celestial

        # Replace NaN values with a very negative number for astrodendro
        clean_data = np.where(np.isfinite(data), data, -1e10)

        # Create dendrogram
        try:
            d = Dendrogram.compute(clean_data, min_value=detection_threshold,
                                 min_npix=5, verbose=False)
        except Exception as e:
            print(f"  Error creating dendrogram: {e}", flush=True)
            return []

        sources = []
        source_name = extract_source_name(filepath)

        print(f"  Found {len(d.leaves)} potential sources", flush=True)

        # Process each leaf (source candidate)
        for i, leaf in enumerate(d.leaves):
            # Get peak value
            peak_value = np.max(clean_data[leaf.indices()])

            # Only keep sources with peak > 10-sigma
            if peak_value >= peak_threshold:
                # Get peak position
                peak_idx = np.unravel_index(np.argmax(clean_data[leaf.indices()]),
                                          clean_data.shape)

                # Convert pixel coordinates to world coordinates
                try:
                    sky_coord = wcs.pixel_to_world(peak_idx[1], peak_idx[0])
                    ra = sky_coord.ra.deg
                    dec = sky_coord.dec.deg
                except Exception as e:
                    print(f"    Warning: Could not convert coordinates for source {i+1}: {e}", flush=True)
                    ra = np.nan
                    dec = np.nan

                # Calculate signal-to-noise ratio
                snr = peak_value / noise_level

                sources.append({
                    'source_field': source_name,
                    'source_number': i + 1,
                    'source_id': f"{source_name}_{i+1:03d}",
                    'ra_deg': ra,
                    'dec_deg': dec,
                    'peak_flux': peak_value,
                    'snr': snr,
                    'noise_level': noise_level,
                    'peak_x_pix': peak_idx[1],
                    'peak_y_pix': peak_idx[0],
                    'npix': len(leaf.indices()[0])
                })

                print(f"    Source {i+1}: SNR={snr:.1f}, Peak={peak_value:.6e}, RA={ra:.6f}, Dec={dec:.6f}", flush=True)

        print(f"  Kept {len(sources)} sources above 10σ threshold", flush=True)
        return sources

def catalog_all_sources():
    """Main function to catalog all sources from FITS files"""

    # Find all FITS files
    fits_pattern = "/orange/adamginsburg/salt/dihca/*quick_cont*fits"
    fits_files = glob.glob(fits_pattern)

    print(f"Found {len(fits_files)} FITS files to process", flush=True)

    all_sources = []
    processed_count = 0
    batch_size = 10  # Process in smaller batches to avoid memory issues

    # Process each FITS file
    for i, fits_file in enumerate(sorted(fits_files)):
        print(f"\nProcessing {i+1}/{len(fits_files)}: {os.path.basename(fits_file)}", flush=True)
        sources = process_fits_file(fits_file)
        all_sources.extend(sources)
        processed_count += 1

        # Save intermediate results every batch_size files
        if processed_count % batch_size == 0:
            print(f"  Processed {processed_count} files so far, found {len(all_sources)} total sources", flush=True)

    if not all_sources:
        print("No sources found!", flush=True)
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
    table.meta['description'] = 'Source catalog from DIHCA continuum images'
    table.meta['detection_threshold'] = '3-sigma'
    table.meta['peak_threshold'] = '10-sigma'
    table.meta['noise_estimator'] = 'mad_std'

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
    # Suppress some warnings
    warnings.filterwarnings('ignore', category=RuntimeWarning)

    catalog = catalog_all_sources()