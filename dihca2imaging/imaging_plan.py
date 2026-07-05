"""
Generate a CASA imaging plan from the source catalog.

This script groups sources spatially and creates tclean imaging commands
for each group, optimizing the field of view based on source distribution.
"""

import numpy as np
import json
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
import os


def read_catalog(catalog_file):
    """Read the source catalog."""
    return Table.read(catalog_file)


def group_sources_by_field(catalog):
    """Group sources by their parent field."""
    fields = {}
    for row in catalog:
        field_name = row['source_field'].decode() if isinstance(row['source_field'], bytes) else row['source_field']
        if field_name not in fields:
            fields[field_name] = []
        fields[field_name].append(row)
    return fields


def cluster_sources(sources, max_separation_arcsec=5.12):
    """
    Cluster sources that are within max_separation_arcsec of each other.
    Uses hierarchical clustering on angular separations.
    """
    if len(sources) == 1:
        return [sources]

    # Extract coordinates
    coords = SkyCoord(
        ra=[s['ra_deg'] for s in sources] * u.deg,
        dec=[s['dec_deg'] for s in sources] * u.deg
    )

    # Compute pairwise separations in arcseconds
    n = len(coords)
    sep_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            sep = coords[i].separation(coords[j]).to(u.arcsec).value
            sep_matrix[i, j] = sep
            sep_matrix[j, i] = sep

    # Convert to condensed distance matrix for hierarchical clustering
    condensed_dist = squareform(sep_matrix)

    # Perform hierarchical clustering
    Z = linkage(condensed_dist, method='complete')

    # Form flat clusters with threshold
    cluster_ids = fcluster(Z, max_separation_arcsec, criterion='distance')

    # Group sources by cluster
    clusters = {}
    for idx, cluster_id in enumerate(cluster_ids):
        if cluster_id not in clusters:
            clusters[cluster_id] = []
        clusters[cluster_id].append(sources[idx])

    return list(clusters.values())


def compute_group_bounds(group):
    """
    Compute the bounding box and optimal image size for a group of sources.
    Returns the center coordinates, image size, and cell size.
    """
    coords = SkyCoord(
        ra=[s['ra_deg'] for s in group] * u.deg,
        dec=[s['dec_deg'] for s in group] * u.deg
    )

    # Compute center
    ra_center = np.mean([s['ra_deg'] for s in group])
    dec_center = np.mean([s['dec_deg'] for s in group])

    # Compute extent
    ra_min = np.min([s['ra_deg'] for s in group])
    ra_max = np.max([s['ra_deg'] for s in group])
    dec_min = np.min([s['dec_deg'] for s in group])
    dec_max = np.max([s['dec_deg'] for s in group])

    # Convert to arcsec (approximate, using dec_center for cos(dec) correction)
    ra_extent_arcsec = (ra_max - ra_min) * 3600 * np.cos(np.radians(dec_center))
    dec_extent_arcsec = (dec_max - dec_min) * 3600

    # Add padding (20% of extent or minimum 1 arcsec)
    padding = max(1.0, 0.2 * max(ra_extent_arcsec, dec_extent_arcsec))
    ra_extent_arcsec += 2 * padding
    dec_extent_arcsec += 2 * padding

    # Determine image size
    cell_size = 0.01  # arcsec per pixel

    # For groups with 5+ sources spread over larger areas, allow up to 15"x15"
    max_size_arcsec = 5.12  # default max
    if len(group) >= 5:
        max_size_arcsec = 15.0

    # Compute required size
    required_size_arcsec = max(ra_extent_arcsec, dec_extent_arcsec)

    # Clamp to allowed range and round up to nearest power of 2 for efficiency
    if required_size_arcsec <= 5.12:
        fov_arcsec = 5.12
        imsize = 512
    elif required_size_arcsec <= 10.24:
        fov_arcsec = 10.24
        imsize = 1024
    else:
        fov_arcsec = 15.36
        imsize = 1536

    return {
        'ra_center': ra_center,
        'dec_center': dec_center,
        'imsize': [imsize, imsize],
        'cell': f'{cell_size}arcsec',
        'fov_arcsec': fov_arcsec,
        'n_sources': len(group)
    }


def find_ms_file(field_name):
    """Find the measurement set file for a given field."""
    # Base directory for MS files
    base_dir = "/orange/adamginsburg/ALMA_IMF/Uploads/dihca2/original_uvdata"

    # Try exact match first
    ms_file = os.path.join(base_dir, f"{field_name}.ms")
    if os.path.exists(ms_file):
        return ms_file

    # Try without _eb suffix
    if '_eb' in field_name:
        field_base = field_name.rsplit('_eb', 1)[0]
        ms_file = os.path.join(base_dir, f"{field_base}.ms")
        if os.path.exists(ms_file):
            return ms_file

    # Return the expected path even if it doesn't exist
    # (may need to be created or the path may be different)
    return os.path.join(base_dir, f"{field_name}.ms")


def get_field_name_for_casa(field_name):
    """Convert field name to CASA field parameter format."""
    # Remove 'G' prefix if present
    if field_name.startswith('G'):
        field_name = field_name[1:]
    # Remove _eb suffix if present
    if '_eb' in field_name:
        field_name = field_name.rsplit('_eb', 1)[0]
    return field_name


def generate_tclean_commands(field_name, groups):
    """
    Generate CASA tclean commands for all groups in a field.
    Returns a dictionary structure similar to dihca2_tclean_commands.json
    """
    commands = {}

    # Find MS file
    vis_file = find_ms_file(field_name)
    casa_field = get_field_name_for_casa(field_name)

    # Spectral windows to image
    spws = [0, 1, 2, 3]

    for group_idx, group in enumerate(groups, start=1):
        group_name = f"{field_name}_group{group_idx}"
        bounds = compute_group_bounds(group)

        # Create phase center string
        phase_center = f"J2000 {bounds['ra_center']}deg {bounds['dec_center']}deg"

        for spw in spws:
            # Create unique key for this group-spw combination
            key = f"{group_name}_spw{spw}"

            # Get rest frequency based on spw (typical ALMA band 6 setup)
            rest_freqs = {
                0: "218.541687GHz",  # H2CO 3(0,3)-2(0,2)
                1: "220.398684GHz",  # 13CO 2-1
                2: "217.104984GHz",  # SiO 5-4
                3: "219.560358GHz"   # C18O 2-1
            }

            commands[key] = {
                "vis": [vis_file],
                "field": casa_field,
                "spw": spw,
                "imagename": f"{group_name}_spw{spw}",
                "phasecenter": phase_center,
                "imsize": bounds['imsize'],
                "cell": bounds['cell'],
                "nchan": 3840,
                "start": 0,
                "width": 1,
                "outframe": "LSRK",
                "veltype": "radio",
                "restfreq": rest_freqs.get(spw, "220.0GHz"),
                "gridder": "standard",
                "specmode": "cube",
                "deconvolver": "hogbom",
                "weighting": "briggs",
                "robust": 0.5,
                "niter": 10000,
                "threshold": "0.1mJy",
                "interactive": False,
                "savemodel": "none",
                "parallel": False,
                "calcres": True,
                "calcpsf": True,
                "metadata": {
                    "n_sources": bounds['n_sources'],
                    "fov_arcsec": bounds['fov_arcsec'],
                    "source_ids": [
                        s['source_id'].decode() if isinstance(s['source_id'], bytes) else s['source_id']
                        for s in group
                    ]
                }
            }

    return commands


def main():
    """Main execution function."""
    print("Reading source catalog...")
    catalog = read_catalog('dihca_source_catalog.fits')
    print(f"Found {len(catalog)} sources")

    print("\nGrouping sources by field...")
    fields = group_sources_by_field(catalog)
    print(f"Found {len(fields)} unique fields")

    print("\nClustering sources and generating imaging commands...")
    all_commands = {}

    for field_name, sources in fields.items():
        print(f"\n  Processing {field_name}: {len(sources)} sources")

        # Cluster sources within this field
        groups = cluster_sources(sources, max_separation_arcsec=5.12)
        print(f"    Created {len(groups)} groups")

        for i, group in enumerate(groups, start=1):
            bounds = compute_group_bounds(group)
            print(f"      Group {i}: {bounds['n_sources']} sources, "
                  f"{bounds['fov_arcsec']:.2f}\" FOV, "
                  f"{bounds['imsize'][0]}x{bounds['imsize'][1]} pixels")

        # Generate tclean commands
        field_commands = generate_tclean_commands(field_name, groups)
        all_commands.update(field_commands)

    # Write output
    output_file = 'dihca2_imaging_plan.json'
    print(f"\nWriting imaging plan to {output_file}...")
    with open(output_file, 'w') as f:
        json.dump(all_commands, f, indent=2)

    print(f"\nGenerated {len(all_commands)} imaging commands")
    print("Done!")

    # Print summary statistics
    total_groups = len(set(cmd.rsplit('_spw', 1)[0] for cmd in all_commands.keys()))
    print(f"\nSummary:")
    print(f"  Total fields: {len(fields)}")
    print(f"  Total groups: {total_groups}")
    print(f"  Total imaging commands: {len(all_commands)}")
    print(f"  Average groups per field: {total_groups / len(fields):.1f}")


if __name__ == "__main__":
    main()
