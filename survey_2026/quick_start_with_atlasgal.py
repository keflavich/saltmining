"""
Quick start: Generate Phase 1 candidate list with ATLASGAL only
Can be updated when other catalogs arrive
"""

import numpy as np
from pathlib import Path
from astropy.table import Table
from astropy.io import fits

DATA_DIR = Path(__file__).parent / 'data'
DUSTMAP_PATH = Path('/orange/adamginsburg/galactic_plane_surveys/dustmaps/edenhofer_mean_and_std_lbd.fits')

def load_atlasgal():
    """Load ATLASGAL catalog"""
    filepath = DATA_DIR / 'atlasgal_clumps.fits'
    if not filepath.exists():
        print(f"ERROR: ATLASGAL catalog not found at {filepath}")
        return None

    print(f"Loading ATLASGAL from {filepath}")
    atlasgal = Table.read(filepath)
    print(f"  Loaded {len(atlasgal)} sources")

    # Check column names
    print(f"  Available columns: {atlasgal.colnames[:15]}")

    # Map columns to standard names
    col_mapping = {}

    # Find distance column
    for col in ['Dist', 'DIST', 'Distance', 'distance_pc']:
        if col in atlasgal.colnames:
            col_mapping['distance_pc'] = col
            break

    # Find luminosity column
    for col in ['Lbol', 'LBOL', 'Lbol_Lsun', 'L_bol_Lsun']:
        if col in atlasgal.colnames:
            col_mapping['L_bol_Lsun'] = col
            break

    # Find Galactic coordinates
    for col in ['GLON', 'Glon', 'GLON_90']:
        if col in atlasgal.colnames:
            col_mapping['l'] = col
            break

    for col in ['GLAT', 'Glat']:
        if col in atlasgal.colnames:
            col_mapping['b'] = col
            break

    # Rename
    for new_name, old_name in col_mapping.items():
        if old_name in atlasgal.colnames and new_name != old_name:
            atlasgal.rename_column(old_name, new_name)

    return atlasgal

def load_dustmap():
    """Load Edenhofer dust map"""
    print(f"\nLoading Edenhofer dust map from {DUSTMAP_PATH}")
    if not DUSTMAP_PATH.exists():
        print(f"  WARNING: Dust map not found")
        return None

    try:
        hdul = fits.open(DUSTMAP_PATH)
        dust_data = hdul[0].data
        dust_header = hdul[0].header
        hdul.close()
        print(f"  Loaded dust map with shape {dust_data.shape}")
        return dust_data, dust_header
    except Exception as e:
        print(f"  ERROR: {e}")
        return None

def filter_catalog(atlasgal, d_min_kpc=0.4, d_max_kpc=2.0, L_min_Lsun=1e4):
    """Filter by distance and luminosity"""
    print(f"\nApplying filters:")
    print(f"  Distance: {d_min_kpc} < d < {d_max_kpc} kpc")
    print(f"  Luminosity: L > {L_min_Lsun:.0e} Lsun")

    # Distance filter
    if 'distance_pc' not in atlasgal.colnames:
        print("  WARNING: No distance column found")
        return atlasgal

    d_kpc = atlasgal['distance_pc'] / 1000.0
    mask = (d_kpc >= d_min_kpc) & (d_kpc <= d_max_kpc)

    # Luminosity filter (if available)
    if 'L_bol_Lsun' in atlasgal.colnames:
        valid_L = ~np.isnan(atlasgal['L_bol_Lsun'])
        mask &= (valid_L & (atlasgal['L_bol_Lsun'] >= L_min_Lsun))
        print(f"  Used luminosity filter")
    else:
        print(f"  WARNING: No luminosity column; filtering by distance only")

    filtered = atlasgal[mask]
    print(f"  Result: {len(filtered)} / {len(atlasgal)} sources")

    return filtered

def print_statistics(catalog):
    """Print basic statistics"""
    print(f"\n{'='*70}")
    print("STATISTICS")
    print(f"{'='*70}")

    if 'distance_pc' in catalog.colnames:
        d_kpc = catalog['distance_pc'] / 1000.0
        print(f"Distance (kpc):")
        print(f"  Min: {d_kpc.min():.2f}")
        print(f"  Max: {d_kpc.max():.2f}")
        print(f"  Mean: {np.nanmean(d_kpc):.2f}")
        print(f"  Median: {np.nanmedian(d_kpc):.2f}")

    if 'L_bol_Lsun' in catalog.colnames:
        L_valid = catalog['L_bol_Lsun'][~np.isnan(catalog['L_bol_Lsun'])]
        if len(L_valid) > 0:
            print(f"Luminosity (Lsun):")
            print(f"  N with L: {len(L_valid)}")
            print(f"  Min: {L_valid.min():.2e}")
            print(f"  Max: {L_valid.max():.2e}")
            print(f"  Median: {np.nanmedian(L_valid):.2e}")

    # Galactic distribution
    if 'l' in catalog.colnames:
        print(f"Galactic coordinates:")
        print(f"  l range: {catalog['l'].min():.1f}° to {catalog['l'].max():.1f}°")
        print(f"  b range: {catalog['b'].min():.1f}° to {catalog['b'].max():.1f}°")

def main():
    print("\n" + "="*70)
    print("QUICK START: PHASE 1 WITH ATLASGAL")
    print("="*70)
    print("This generates an initial Phase 1 catalog with ATLASGAL data only.")
    print("Will be updated when other catalogs (RMS, SEDIGISM, WGO, Hi-GAL) download.")

    # Load
    atlasgal = load_atlasgal()
    if atlasgal is None:
        print("\nERROR: Cannot proceed without ATLASGAL")
        return

    dust_data = load_dustmap()

    # Filter
    filtered = filter_catalog(atlasgal)

    # Statistics
    print_statistics(filtered)

    # Save
    output_path = DATA_DIR / 'myso_candidates_phase1_atlasgal_only.fits'
    filtered.write(output_path, overwrite=True)
    print(f"\n✓ Saved to {output_path}")

    print(f"\n{'='*70}")
    print("NEXT STEPS")
    print(f"{'='*70}")
    print("""
1. Check download status: ls -lh data/
2. Once other catalogs download, run: python load_and_crossmatch.py
   This will merge all catalogs and produce the full Phase 1 table
3. For any missing catalogs, download manually from:
   - RMS: http://rms.leeds.ac.uk/
   - SEDIGISM: http://sedigism.mpifr-bonn.mpg.de/
   - WISE Green Objects: https://ui.adsabs.harvard.edu/abs/2023yCat..22640024Z/abstract
   - Hi-GAL: https://ui.adsabs.harvard.edu/abs/2013yCat..35510100B/abstract
""")

if __name__ == '__main__':
    main()
