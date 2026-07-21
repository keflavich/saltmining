#!/usr/bin/env python3
"""
Script to inspect MS files and extract metadata for imaging parameters.
This script examines all MS files in the DIHCA2 dataset to extract:
- Spectral window information
- Field names and coordinates
- UV range (for determining image size)
- Channel information
- Other relevant metadata for tclean parameters
"""

import os
import sys
import json
import glob
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import argparse


def inspect_ms_with_casa(ms_path):
    """
    Inspect an MS file using CASA tools to extract metadata.
    This function should be run within a CASA environment.
    """
    # Import CASA tools
    from casatools import msmetadata, table, measures

    msmd = msmetadata()
    tb = table()
    me = measures()

    print(f"Inspecting {ms_path}")

    # Open MS metadata
    msmd.open(ms_path)

    # Get basic information
    field_names = msmd.fieldnames()
    field_ids = msmd.fieldsforname(field_names[0])  # Assuming single field per MS
    spw_ids = msmd.spwsforfield(field_ids[0])

    metadata = {
        'ms_path': ms_path,
        'field_names': field_names,
        'field_ids': field_ids.tolist(),
        'spw_ids': spw_ids.tolist(),
        'spectral_windows': {}
    }

    # Get field coordinates
    field_directions = []
    for field_id in field_ids:
        direction = msmd.phasecenter(field_id)
        ra_rad = direction['m0']['value']
        dec_rad = direction['m1']['value']

        # Convert to degrees
        ra_deg = np.degrees(ra_rad)
        dec_deg = np.degrees(dec_rad)

        field_directions.append({
            'field_id': int(field_id),
            'ra_deg': ra_deg,
            'dec_deg': dec_deg,
            'ra_rad': ra_rad,
            'dec_rad': dec_rad
        })

    metadata['field_directions'] = field_directions

    # Get spectral window information
    for spw_id in spw_ids:
        spw_info = {
            'spw_id': int(spw_id),
            'nchan': int(msmd.nchan(spw_id)),
            'ref_freq': msmd.reffreq(spw_id)['m0']['value'],  # Hz
            'chan_freq': msmd.chanfreqs(spw_id).tolist(),  # Hz
            'chan_width': msmd.chanwidths(spw_id).tolist(),  # Hz
            'bandwidth': msmd.bandwidths(spw_id),  # Hz
            'name': msmd.namesforspws(spw_id)[0] if msmd.namesforspws(spw_id) else f'spw{spw_id}'
        }

        # Calculate frequency range
        frequencies = np.array(spw_info['chan_freq'])
        spw_info['freq_min'] = float(frequencies.min())
        spw_info['freq_max'] = float(frequencies.max())
        spw_info['freq_center'] = float(np.mean(frequencies))

        metadata['spectral_windows'][str(spw_id)] = spw_info

    # Get UV range information
    tb.open(ms_path)
    uvw = tb.getcol('UVW')
    tb.close()

    # Calculate UV distances
    uv_dist = np.sqrt(uvw[0]**2 + uvw[1]**2)  # meters

    # Convert to lambda at reference frequency (use first spw)
    ref_freq = metadata['spectral_windows'][str(spw_ids[0])]['ref_freq']  # Hz
    wavelength = 2.998e8 / ref_freq  # meters
    uv_dist_lambda = uv_dist / wavelength

    metadata['uv_range'] = {
        'min_uv_m': float(uv_dist.min()),
        'max_uv_m': float(uv_dist.max()),
        'min_uv_lambda': float(uv_dist_lambda.min()),
        'max_uv_lambda': float(uv_dist_lambda.max()),
        'ref_freq_hz': ref_freq,
        'wavelength_m': wavelength
    }

    # Calculate suggested image parameters
    # Image size should be ~3x the primary beam FWHM
    # For ALMA: theta_PB ≈ 1.13 * lambda/D where D is dish diameter (12m for ALMA)
    dish_diameter = 12.0  # meters
    primary_beam_fwhm_rad = 1.13 * wavelength / dish_diameter
    primary_beam_fwhm_arcsec = np.degrees(primary_beam_fwhm_rad) * 3600

    # Suggested image size (3x primary beam)
    suggested_imsize_arcsec = 3 * primary_beam_fwhm_arcsec

    # Cell size should be ~lambda/max_baseline/5
    max_baseline = uv_dist.max()
    suggested_cell_arcsec = np.degrees(wavelength / max_baseline / 5) * 3600

    # Round to reasonable values
    suggested_imsize_pixels = int(np.ceil(suggested_imsize_arcsec / suggested_cell_arcsec))
    # Round to power of 2 for efficiency
    suggested_imsize_pixels = int(2**np.ceil(np.log2(suggested_imsize_pixels)))

    # Apply maximum image size constraint
    max_imsize = 4096
    if suggested_imsize_pixels > max_imsize:
        # Calculate how much we need to scale up the cell size
        scale_factor = suggested_imsize_pixels / max_imsize
        # Limit the scale factor to 3x maximum
        scale_factor = min(scale_factor, 3.0)

        # Apply the scaling
        suggested_cell_arcsec *= scale_factor
        suggested_imsize_pixels = max_imsize

        print(f"  Warning: Image size reduced from {int(np.ceil(suggested_imsize_arcsec / (suggested_cell_arcsec / scale_factor)))} to {max_imsize} pixels")
        print(f"  Cell size increased by {scale_factor:.2f}x to {suggested_cell_arcsec:.3f} arcsec")

    metadata['suggested_imaging_params'] = {
        'primary_beam_fwhm_arcsec': primary_beam_fwhm_arcsec,
        'suggested_imsize_arcsec': suggested_imsize_arcsec,
        'suggested_cell_arcsec': suggested_cell_arcsec,
        'suggested_imsize_pixels': suggested_imsize_pixels,
        'max_baseline_m': float(max_baseline)
    }

    msmd.close()

    return metadata


def inspect_all_ms_files(data_dir, output_file=None, update_existing=True, force_update=False, limit=None, limit_incomplete=None):
    """
    Inspect all MS files in the data directory and collect metadata.
    Can work incrementally by updating existing metadata file.
    """
    ms_files = glob.glob(os.path.join(data_dir, "*.ms"))
    ms_files.sort()

    if limit:
        ms_files = ms_files[:limit]
        print(f"Limited to first {limit} MS files")

    print(f"Found {len(ms_files)} MS files")

    # Load existing metadata if available and requested
    all_metadata = {}
    if update_existing and output_file and os.path.exists(output_file):
        try:
            with open(output_file, 'r') as f:
                all_metadata = json.load(f)
            print(f"Loaded existing metadata with {len(all_metadata)} entries")
        except Exception as e:
            print(f"Could not load existing metadata: {e}")
            all_metadata = {}

    processed = 0
    skipped = 0

    for ms_file in ms_files:
        ms_name = os.path.basename(ms_file)

        # Skip if already processed (unless forcing update)
        if ms_name in all_metadata and not force_update:
            # Check if this entry has the new image size constraints
            if ('suggested_imaging_params' in all_metadata[ms_name] and
                all_metadata[ms_name]['suggested_imaging_params'].get('suggested_imsize_pixels', 0) <= 4096):
                print(f"Skipping {ms_name} (already processed with size constraints)")
                skipped += 1
                continue
            else:
                print(f"Updating {ms_name} (needs size constraints applied)")
        else:
            status = "forced update" if force_update and ms_name in all_metadata else "new entry"
            print(f"Processing {ms_name} ({status})")

        try:
            metadata = inspect_ms_with_casa(ms_file)
            all_metadata[ms_name] = metadata
            processed += 1

            # Save incrementally every 5 files to avoid losing work
            if output_file and processed % 5 == 0:
                with open(output_file, 'w') as f:
                    json.dump(all_metadata, f, indent=2)
                print(f"  Saved incremental progress ({processed} processed, {skipped} skipped)")

        except Exception as e:
            print(f"Error processing {ms_name}: {e}")
            continue

    print(f"\nCompleted: {processed} processed, {skipped} skipped, {len(all_metadata)} total entries")
    return all_metadata


def generate_tclean_parameters(metadata_dict, output_dir="/orange/adamginsburg/salt/dihca2imaging"):
    """
    Generate tclean parameters based on the MS metadata, following ACES pipeline format.
    """
    tclean_commands = {}

    for ms_name, metadata in metadata_dict.items():
        field_name = metadata['field_names'][0] if metadata['field_names'] else ms_name.replace('.ms', '')

        # Create entry for this MS
        ms_key = ms_name.replace('.ms', '')
        tclean_commands[ms_key] = {
            'tclean_cube_pars': {},
            'field_name': field_name,
            'ms_path': metadata['ms_path']
        }

        # Generate parameters for each spectral window
        for spw_id_str, spw_info in metadata['spectral_windows'].items():
            spw_id = int(spw_id_str)

            # Base tclean parameters
            base_params = {
                'vis': [metadata['ms_path']],
                'field': field_name,
                'spw': spw_id,
                'imagename': f"{field_name}_spw{spw_id}",
                'imsize': [metadata['suggested_imaging_params']['suggested_imsize_pixels']] * 2,
                'cell': f"{metadata['suggested_imaging_params']['suggested_cell_arcsec']:.3f}arcsec",
                'nchan': spw_info['nchan'],
                'start': 0,
                'width': 1,
                'outframe': 'LSRK',
                'veltype': 'radio',
                'restfreq': f"{spw_info['ref_freq']/1e9:.6f}GHz",
                'gridder': 'mosaic',
                'specmode': 'cube',
                'deconvolver': 'hogbom',
                'weighting': 'briggs',
                'robust': 0.5,
                'niter': 10000,
                'threshold': '0.1mJy',
                'interactive': False,
                'savemodel': 'none',
                'parallel': False,
                'calcres': True,
                'calcpsf': True
            }

            # Store parameters
            tclean_commands[ms_key]['tclean_cube_pars'][spw_id_str] = base_params

    return tclean_commands


def main():
    parser = argparse.ArgumentParser(description='Inspect MS files and generate imaging parameters')
    parser.add_argument('--data-dir', default='/orange/adamginsburg/ALMA_IMF/Uploads/dihca2/original_uvdata',
                        help='Directory containing MS files')
    parser.add_argument('--output-dir', default='/orange/adamginsburg/salt/dihca2imaging',
                        help='Output directory for JSON files')
    parser.add_argument('--output-file', default='dihca2_metadata.json',
                        help='Output JSON file name')
    parser.add_argument('--force-update', action='store_true',
                        help='Force update of all entries, even if already processed')
    parser.add_argument('--limit', type=int,
                        help='Limit processing to first N MS files (for testing)')
    parser.add_argument('--limit-incomplete', type=int,
                        help='Limit processing to next N incomplete MS files (missing or needing size constraints)')

    args = parser.parse_args()

    # Check if we're running in CASA
    from casatools import msmetadata

    # Inspect all MS files
    print(f"Inspecting MS files in {args.data_dir}")
    metadata_file = os.path.join(args.output_dir, args.output_file)
    metadata = inspect_all_ms_files(args.data_dir,
                                   output_file=metadata_file,
                                   update_existing=True,
                                   force_update=args.force_update,
                                   limit=args.limit,
                                   limit_incomplete=args.limit_incomplete)

    # Generate tclean parameters
    print("Generating tclean parameters...")
    tclean_params = generate_tclean_parameters(metadata)

    # Save metadata
    os.makedirs(args.output_dir, exist_ok=True)
    metadata_file = os.path.join(args.output_dir, args.output_file)

    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=2)

    print(f"Saved metadata to {metadata_file}")

    # Save tclean parameters
    tclean_file = os.path.join(args.output_dir, 'dihca2_tclean_commands.json')
    with open(tclean_file, 'w') as f:
        json.dump(tclean_params, f, indent=2)

    print(f"Saved tclean parameters to {tclean_file}")

    # Print summary
    print(f"\nSummary:")
    print(f"  - Processed {len(metadata)} MS files")
    print(f"  - Total spectral windows: {sum(len(ms['spectral_windows']) for ms in metadata.values())}")
    print(f"  - Field names: {set(ms['field_names'][0] for ms in metadata.values() if ms['field_names'])}")


if __name__ == '__main__':
    main()
