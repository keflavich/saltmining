#!/usr/bin/env python
"""
Utility script to query and analyze the DIHCA2 imaging plan.

Usage examples:
    # Show all groups for a specific field
    python query_imaging_plan.py --field G008.672-00.682

    # Show groups with large number of sources
    python query_imaging_plan.py --min-sources 10

    # Show groups with specific image size
    python query_imaging_plan.py --imsize 1024

    # Export specific field commands to separate JSON
    python query_imaging_plan.py --field G345.493+01.469 --export
"""

import json
import argparse
import sys
from collections import defaultdict


def load_imaging_plan(filename='dihca2_imaging_plan.json'):
    """Load the imaging plan JSON file."""
    with open(filename, 'r') as f:
        return json.load(f)


def get_unique_groups(imaging_plan):
    """Extract unique group information from the imaging plan."""
    groups = {}
    for key, cmd in imaging_plan.items():
        # Extract group name (remove _spwN suffix)
        group_key = key.rsplit('_spw', 1)[0]
        if group_key not in groups:
            groups[group_key] = {
                'field': cmd['field'],
                'n_sources': cmd['metadata']['n_sources'],
                'fov_arcsec': cmd['metadata']['fov_arcsec'],
                'imsize': cmd['imsize'][0],
                'source_ids': cmd['metadata']['source_ids'],
                'phasecenter': cmd['phasecenter'],
                'ra_center': float(cmd['phasecenter'].split()[1].rstrip('deg')),
                'dec_center': float(cmd['phasecenter'].split()[2].rstrip('deg'))
            }
    return groups


def query_by_field(groups, field_name):
    """Find all groups for a specific field."""
    matching = []
    for group_key, group_data in groups.items():
        # Match field name (with or without 'G' prefix)
        field = group_data['field']
        if field_name.lower() in field.lower() or field_name.lower() in group_key.lower():
            matching.append((group_key, group_data))
    return matching


def query_by_min_sources(groups, min_sources):
    """Find all groups with at least min_sources sources."""
    matching = []
    for group_key, group_data in groups.items():
        if group_data['n_sources'] >= min_sources:
            matching.append((group_key, group_data))
    return sorted(matching, key=lambda x: x[1]['n_sources'], reverse=True)


def query_by_imsize(groups, imsize):
    """Find all groups with specific image size."""
    matching = []
    for group_key, group_data in groups.items():
        if group_data['imsize'] == imsize:
            matching.append((group_key, group_data))
    return matching


def print_group_info(group_key, group_data, verbose=False):
    """Print information about a group."""
    print(f"\n{group_key}")
    print(f"  Field: {group_data['field']}")
    print(f"  Sources: {group_data['n_sources']}")
    print(f"  FOV: {group_data['fov_arcsec']:.2f}\"")
    print(f"  Image size: {group_data['imsize']}×{group_data['imsize']} pixels")
    print(f"  Phase center: {group_data['phasecenter']}")

    if verbose:
        print(f"  Source IDs:")
        for src_id in group_data['source_ids']:
            print(f"    - {src_id}")


def export_field_commands(imaging_plan, field_name, output_file=None):
    """Export all commands for a specific field to a separate JSON file."""
    matching_commands = {}
    for key, cmd in imaging_plan.items():
        field = cmd['field']
        if field_name.lower() in field.lower() or field_name.lower() in key.lower():
            matching_commands[key] = cmd

    if not output_file:
        # Generate output filename from field name
        clean_field = field_name.replace('G', '').replace('+', 'p').replace('-', 'm').replace('.', '_')
        output_file = f'imaging_plan_{clean_field}.json'

    with open(output_file, 'w') as f:
        json.dump(matching_commands, f, indent=2)

    print(f"\nExported {len(matching_commands)} commands to {output_file}")
    return output_file


def show_statistics(groups):
    """Show overall statistics about the imaging plan."""
    import numpy as np

    n_sources = [g['n_sources'] for g in groups.values()]
    imsizes = [g['imsize'] for g in groups.values()]
    fovs = [g['fov_arcsec'] for g in groups.values()]

    # Count by field
    fields = defaultdict(int)
    for g in groups.values():
        fields[g['field']] += 1

    print("\n" + "="*60)
    print("IMAGING PLAN STATISTICS")
    print("="*60)
    print(f"Total groups: {len(groups)}")
    print(f"Total fields: {len(fields)}")
    print()

    print("Sources per group:")
    print(f"  Min: {min(n_sources)}")
    print(f"  Max: {max(n_sources)}")
    print(f"  Mean: {np.mean(n_sources):.1f}")
    print(f"  Median: {np.median(n_sources):.0f}")
    print()

    print("Image size distribution:")
    size_counts = defaultdict(int)
    for size in imsizes:
        size_counts[size] += 1
    for size in sorted(size_counts.keys()):
        fov = size * 0.01
        pct = 100 * size_counts[size] / len(groups)
        print(f"  {size}×{size} pixels ({fov:.2f}\" FOV): {size_counts[size]} groups ({pct:.1f}%)")
    print()

    print("Top 10 fields by number of groups:")
    sorted_fields = sorted(fields.items(), key=lambda x: x[1], reverse=True)[:10]
    for field, count in sorted_fields:
        print(f"  {field}: {count} groups")


def main():
    parser = argparse.ArgumentParser(
        description='Query and analyze the DIHCA2 imaging plan',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument('--field', type=str, help='Query by field name')
    parser.add_argument('--min-sources', type=int, help='Query groups with at least N sources')
    parser.add_argument('--imsize', type=int, help='Query groups with specific image size')
    parser.add_argument('--export', action='store_true', help='Export field commands to separate JSON')
    parser.add_argument('--verbose', '-v', action='store_true', help='Show detailed information')
    parser.add_argument('--stats', action='store_true', help='Show overall statistics')
    parser.add_argument('--input', type=str, default='dihca2_imaging_plan.json',
                        help='Input imaging plan JSON file')

    args = parser.parse_args()

    # Load imaging plan
    try:
        imaging_plan = load_imaging_plan(args.input)
    except FileNotFoundError:
        print(f"Error: Could not find {args.input}")
        sys.exit(1)

    # Get unique groups
    groups = get_unique_groups(imaging_plan)

    # Handle different query modes
    if args.stats:
        show_statistics(groups)
        return

    matching_groups = None

    if args.field:
        matching_groups = query_by_field(groups, args.field)
        if args.export:
            export_field_commands(imaging_plan, args.field)
        print(f"\nFound {len(matching_groups)} groups for field matching '{args.field}':")

    elif args.min_sources:
        matching_groups = query_by_min_sources(groups, args.min_sources)
        print(f"\nFound {len(matching_groups)} groups with {args.min_sources}+ sources:")

    elif args.imsize:
        matching_groups = query_by_imsize(groups, args.imsize)
        print(f"\nFound {len(matching_groups)} groups with image size {args.imsize}×{args.imsize}:")

    else:
        # Default: show statistics
        show_statistics(groups)
        return

    # Print matching groups
    if matching_groups:
        for group_key, group_data in matching_groups:
            print_group_info(group_key, group_data, verbose=args.verbose)


if __name__ == '__main__':
    main()
