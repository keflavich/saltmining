#!/usr/bin/env python
"""
Export imaging groups as DS9 region files.

This script creates DS9-compatible region files showing:
- Group coverage areas (boxes representing the FOV)
- Source positions
- Group labels

Usage:
    # Export all groups to one region file
    python export_regions.py

    # Export regions for a specific field
    python export_regions.py --field G008.672-00.682

    # Export with source positions
    python export_regions.py --show-sources

    # Export separate file per field
    python export_regions.py --separate-fields
"""

import json
import argparse
import sys
from collections import defaultdict
from astropy.table import Table
import numpy as np


def load_imaging_plan(filename='dihca2_imaging_plan.json'):
    """Load the imaging plan JSON file."""
    with open(filename, 'r') as f:
        return json.load(f)


def load_source_catalog(filename='dihca_source_catalog.fits'):
    """Load the source catalog."""
    try:
        return Table.read(filename)
    except FileNotFoundError:
        print(f"Warning: Source catalog {filename} not found. Source positions will not be included.")
        return None


def get_unique_groups(imaging_plan):
    """Extract unique group information from the imaging plan."""
    groups = {}
    for key, cmd in imaging_plan.items():
        # Extract group name (remove _spwN suffix)
        group_key = key.rsplit('_spw', 1)[0]
        if group_key not in groups:
            # Parse phase center
            phasecenter = cmd['phasecenter']
            parts = phasecenter.split()
            ra_deg = float(parts[1].rstrip('deg'))
            dec_deg = float(parts[2].rstrip('deg'))

            groups[group_key] = {
                'field': cmd['field'],
                'n_sources': cmd['metadata']['n_sources'],
                'fov_arcsec': cmd['metadata']['fov_arcsec'],
                'imsize': cmd['imsize'][0],
                'source_ids': cmd['metadata']['source_ids'],
                'ra_deg': ra_deg,
                'dec_deg': dec_deg
            }
    return groups


def write_ds9_region_file(groups, filename, source_catalog=None, field_filter=None):
    """
    Write groups to a DS9 region file.

    Parameters:
    -----------
    groups : dict
        Dictionary of group information
    filename : str
        Output filename
    source_catalog : Table, optional
        Source catalog for plotting individual sources
    field_filter : str, optional
        Only include groups from this field
    """
    with open(filename, 'w') as f:
        # Header
        f.write("# Region file format: DS9 version 4.1\n")
        f.write("global color=green dashlist=8 3 width=2 font=\"helvetica 10 normal roman\" ")
        f.write("select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
        f.write("fk5\n")

        # Filter groups if requested
        if field_filter:
            filtered_groups = {k: v for k, v in groups.items()
                             if field_filter.lower() in k.lower() or
                                field_filter.lower() in v['field'].lower()}
        else:
            filtered_groups = groups

        # Color palette for different fields
        colors = ['green', 'cyan', 'magenta', 'yellow', 'red', 'blue', 'white']
        field_colors = {}
        color_idx = 0

        # Assign colors to fields
        for group_data in filtered_groups.values():
            field = group_data['field']
            if field not in field_colors:
                field_colors[field] = colors[color_idx % len(colors)]
                color_idx += 1

        # Write group boxes
        for group_key, group_data in filtered_groups.items():
            ra = group_data['ra_deg']
            dec = group_data['dec_deg']
            fov_arcsec = group_data['fov_arcsec']
            n_sources = group_data['n_sources']
            field = group_data['field']
            color = field_colors[field]

            # Convert FOV to degrees
            fov_deg = fov_arcsec / 3600.0

            # Write box region (RA, Dec, width, height in degrees)
            # Box is centered at phase center with size matching FOV
            f.write(f"box({ra:.8f},{dec:.8f},{fov_deg:.8f}d,{fov_deg:.8f}d,0) # ")
            f.write(f"color={color} width=2 ")
            f.write(f'text="{group_key} ({n_sources} src)"\n')

        # Write source positions if catalog provided
        if source_catalog is not None:
            f.write("\n# Individual sources\n")

            # Create a set of source IDs that are in our groups
            included_sources = set()
            for group_data in filtered_groups.values():
                included_sources.update(group_data['source_ids'])

            # Plot sources
            for row in source_catalog:
                source_id = row['source_id'].decode() if isinstance(row['source_id'], bytes) else row['source_id']

                if source_id in included_sources:
                    ra = row['ra_deg']
                    dec = row['dec_deg']

                    # Find which group this source belongs to
                    source_color = 'white'
                    for group_key, group_data in filtered_groups.items():
                        if source_id in group_data['source_ids']:
                            source_color = field_colors[group_data['field']]
                            break

                    # Write point region
                    f.write(f"point({ra:.8f},{dec:.8f}) # point=x color={source_color} ")
                    f.write(f'text="{source_id}"\n')

    print(f"Wrote {len(filtered_groups)} groups to {filename}")


def write_separate_field_regions(groups, source_catalog=None, output_dir='.'):
    """Write separate region files for each field."""
    import os

    # Group by field
    fields = defaultdict(dict)
    for group_key, group_data in groups.items():
        field = group_data['field']
        fields[field][group_key] = group_data

    # Write one file per field
    for field, field_groups in fields.items():
        # Clean field name for filename
        clean_field = field.replace('+', 'p').replace('-', 'm').replace('.', '_')
        filename = os.path.join(output_dir, f'imaging_groups_{clean_field}.reg')

        write_ds9_region_file(field_groups, filename, source_catalog=source_catalog)

    print(f"\nWrote {len(fields)} field region files to {output_dir}")


def write_crtf_region_file(groups, filename, field_filter=None):
    """
    Write groups to a CASA CRTF region file.

    Parameters:
    -----------
    groups : dict
        Dictionary of group information
    filename : str
        Output filename
    field_filter : str, optional
        Only include groups from this field
    """
    with open(filename, 'w') as f:
        # Header
        f.write("#CRTFv0 CASA Region Text Format version 0\n")
        f.write("# Imaging groups for DIHCA2\n")
        f.write("global coord=J2000\n\n")

        # Filter groups if requested
        if field_filter:
            filtered_groups = {k: v for k, v in groups.items()
                             if field_filter.lower() in k.lower() or
                                field_filter.lower() in v['field'].lower()}
        else:
            filtered_groups = groups

        # Write group boxes
        for group_key, group_data in filtered_groups.items():
            ra = group_data['ra_deg']
            dec = group_data['dec_deg']
            fov_arcsec = group_data['fov_arcsec']
            n_sources = group_data['n_sources']

            # Convert coordinates to HMS/DMS format for CASA
            from astropy.coordinates import SkyCoord
            from astropy import units as u

            coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
            ra_str = coord.ra.to_string(unit=u.hourangle, sep=':', precision=3, pad=True)
            dec_str = coord.dec.to_string(unit=u.deg, sep='.', precision=2, alwayssign=True, pad=True)

            # Write box region
            half_fov = fov_arcsec / 2.0
            f.write(f"box[[{ra_str}, {dec_str}], [{fov_arcsec}arcsec, {fov_arcsec}arcsec]] ")
            f.write(f"# {group_key} ({n_sources} sources)\n")

    print(f"Wrote {len(filtered_groups)} groups to {filename} (CASA CRTF format)")


def main():
    parser = argparse.ArgumentParser(
        description='Export imaging groups as region files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--field', type=str,
                        help='Export only groups from this field')
    parser.add_argument('--show-sources', action='store_true',
                        help='Include individual source positions')
    parser.add_argument('--separate-fields', action='store_true',
                        help='Write separate region file for each field')
    parser.add_argument('--format', choices=['ds9', 'crtf', 'both'], default='ds9',
                        help='Output format (DS9 or CASA CRTF)')
    parser.add_argument('--input', type=str, default='dihca2_imaging_plan.json',
                        help='Input imaging plan JSON file')
    parser.add_argument('--catalog', type=str, default='dihca_source_catalog.fits',
                        help='Source catalog FITS file')
    parser.add_argument('--output', type=str, default='dihca2_imaging_groups.reg',
                        help='Output region filename (ignored with --separate-fields)')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Output directory for separate field files')

    args = parser.parse_args()

    # Load imaging plan
    print(f"Loading imaging plan from {args.input}...")
    try:
        imaging_plan = load_imaging_plan(args.input)
    except FileNotFoundError:
        print(f"Error: Could not find {args.input}")
        sys.exit(1)

    # Get unique groups
    groups = get_unique_groups(imaging_plan)
    print(f"Found {len(groups)} unique groups")

    # Load source catalog if needed
    source_catalog = None
    if args.show_sources:
        print(f"Loading source catalog from {args.catalog}...")
        source_catalog = load_source_catalog(args.catalog)

    # Export regions
    if args.separate_fields:
        # Write separate files
        if args.format in ['ds9', 'both']:
            write_separate_field_regions(groups, source_catalog=source_catalog,
                                        output_dir=args.output_dir)
    else:
        # Write single file
        if args.format in ['ds9', 'both']:
            write_ds9_region_file(groups, args.output,
                                 source_catalog=source_catalog,
                                 field_filter=args.field)

        if args.format in ['crtf', 'both']:
            crtf_file = args.output.replace('.reg', '.crtf')
            write_crtf_region_file(groups, crtf_file, field_filter=args.field)

    print("\nDone!")


if __name__ == '__main__':
    main()
