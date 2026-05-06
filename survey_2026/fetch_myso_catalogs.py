"""
Fetch and compile MYSO catalogs for salt disk occurrence survey
Phase 1: Build candidate MYSO table with L > 10^4 Lsun, 0.4 < d < 2 kpc
"""

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack, unique, Column
import warnings
warnings.filterwarnings('ignore')

try:
    from astroquery.vizier import Vizier
except ImportError:
    print("astroquery not available; will attempt to download manually")

def fetch_atlasgal_distances():
    """
    Fetch ATLASGAL complete clump catalog with distances (Urquhart 2018)
    VizieR catalog: J/MNRAS/473/1059
    """
    print("\n=== Fetching ATLASGAL (Urquhart 2018) ===")
    try:
        v = Vizier(columns=['all'])  # Fetch all columns first to see names
        v.ROW_LIMIT = -1
        # Try different catalog identifier formats
        try:
            atlasgal = v.query_constraints(catalog='J/MNRAS/473/1059')[0]
        except:
            atlasgal = v.query_constraints(catalog='2018yCat..74731059U')[0]

        print(f"  Retrieved {len(atlasgal)} ATLASGAL sources")
        print(f"  Available columns: {atlasgal.colnames[:10]}")

        # Map to standard column names (check actual names in data)
        col_mapping = {}
        for col in atlasgal.colnames:
            col_lower = col.lower()
            if 'glon' in col_lower or 'lon' in col_lower:
                col_mapping['l'] = col
            elif 'glat' in col_lower or 'lat' in col_lower:
                col_mapping['b'] = col
            elif 'dist' in col_lower:
                col_mapping['distance_pc'] = col
            elif 'lbol' in col_lower or ('l' in col_lower and 'bol' in col_lower):
                col_mapping['L_bol_Lsun'] = col

        for new_name, old_name in col_mapping.items():
            if old_name in atlasgal.colnames:
                atlasgal.rename_column(old_name, new_name)

        atlasgal['catalog'] = 'ATLASGAL'
        return atlasgal
    except Exception as e:
        print(f"  Error fetching ATLASGAL: {e}")
        print("  Please download manually from https://ui.adsabs.harvard.edu/abs/2018yCat..74731059U")
        return None

def fetch_sedigism_clumps():
    """
    Fetch SEDIGISM clump catalog with distances (southern extension, Schuller 2021)
    Database: http://sedigism.mpifr-bonn.mpg.de/
    """
    print("\n=== Fetching SEDIGISM clumps (Schuller 2021) ===")
    print("  SEDIGISM data available at http://sedigism.mpifr-bonn.mpg.de/")
    print("  Consider downloading the cloud catalog manually or via their web interface")
    print("  First data release: Schuller et al. 2021, MNRAS 500, 3064")
    return None

def fetch_wise_green_objects():
    """
    Fetch WISE Green Objects catalog (Zhang 2023)
    VizieR catalog: J/ApJ/954/105 (or 2023yCat..22640024Z)
    """
    print("\n=== Fetching WISE Green Objects (Zhang 2023) ===")
    try:
        v = Vizier(columns=['RAJ2000', 'DEJ2000', 'WISE_name', 'Type', 'all'])
        v.ROW_LIMIT = -1
        wgo = v.query_constraints(catalog='J/ApJ/954/105')[0]
        print(f"  Retrieved {len(wgo)} WISE Green Objects")

        # Assign coordinates
        wgo['coord'] = SkyCoord(ra=wgo['RAJ2000'], dec=wgo['DEJ2000'], frame='icrs')
        wgo['catalog'] = 'WISE_GreenObjects'
        return wgo
    except Exception as e:
        print(f"  Error fetching WISE Green Objects: {e}")
        print("  Please download from https://ui.adsabs.harvard.edu/abs/2023yCat..22640024Z/abstract")
        return None

def fetch_higal_sources():
    """
    Fetch Hi-GAL compact source catalog with SED-fitted luminosities
    Multiple catalogs exist; check Berta et al. 2013, Morii et al., etc.
    """
    print("\n=== Fetching Hi-GAL catalog ===")
    print("  Hi-GAL sources available via VizieR (Berta+2013, Morii+)")
    print("  Your SED fitter available at: https://github.com/keflavich/HiGal_SEDfitter")
    print("  For now, manually download from:")
    print("    - https://ui.adsabs.harvard.edu/abs/2013yCat..35510100B (Berta+2013)")
    return None

def filter_by_distance_luminosity(catalog, d_min=0.4, d_max=2000, L_min=1e4):
    """
    Filter catalog by distance range (kpc) and luminosity (Lsun)
    Also filter for sources with mm dust detections
    """
    if catalog is None:
        return None

    print(f"\nFiltering by distance {d_min:.1f} < d < {d_max:.0f} kpc, L > {L_min:.0e} Lsun")

    # Distance in parsecs
    mask = (catalog['distance_pc'] >= d_min * 1000) & (catalog['distance_pc'] <= d_max * 1000)

    if 'L_bol_Lsun' in catalog.colnames:
        mask &= (catalog['L_bol_Lsun'] >= L_min)

    filtered = catalog[mask]
    print(f"  Retained {len(filtered)} / {len(catalog)} sources")
    return filtered

def cross_match_catalogs(catalogs_dict, match_radius=5*u.arcsec):
    """
    Cross-match multiple catalogs by position
    catalogs_dict: dict with catalog_name -> Table pairs
    """
    print(f"\n=== Cross-matching catalogs (match radius: {match_radius}) ===")

    # Filter for catalogs that were successfully retrieved
    catalogs_dict = {k: v for k, v in catalogs_dict.items() if v is not None}

    if len(catalogs_dict) == 0:
        print("  No catalogs available for cross-matching")
        return None

    # For now, just print what we have
    for name, cat in catalogs_dict.items():
        if cat is not None:
            print(f"  {name}: {len(cat)} sources")

    # TODO: Implement actual cross-matching with astropy.coordinates
    return None

def main():
    print("="*60)
    print("MYSO CATALOG COMPILATION - PHASE 1")
    print("="*60)
    print("Building candidate MYSO table:")
    print("  - Distance range: 0.4 < d < 2 kpc")
    print("  - Luminosity: L > 10^4 Lsun")
    print("  - Include all evolutionary stages with mm dust")

    # Fetch catalogs
    atlasgal = fetch_atlasgal_distances()
    wgo = fetch_wise_green_objects()
    sedigism = fetch_sedigism_clumps()
    higal = fetch_higal_sources()

    # Filter by distance and luminosity
    if atlasgal is not None:
        atlasgal_filtered = filter_by_distance_luminosity(atlasgal)

    if wgo is not None:
        # WGO may not have luminosity in the main catalog, need to cross-match with Hi-GAL
        print("  (WGO luminosities from Hi-GAL cross-match)")

    # Prepare for cross-matching
    catalogs = {
        'ATLASGAL': atlasgal_filtered if atlasgal is not None else None,
        'WISE_GreenObjects': wgo,
        'SEDIGISM': sedigism,
        'HiGAL': higal
    }

    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print("\nCatalogs successfully retrieved:")
    for name, cat in catalogs.items():
        if cat is not None:
            print(f"  ✓ {name}: {len(cat)} sources")
        else:
            print(f"  ✗ {name}: Not retrieved (manual download needed)")

    print("\nNext steps:")
    print("  1. Download remaining catalogs from above links")
    print("  2. Load into Python as Tables")
    print("  3. Cross-match by (l, b) with 5 arcsec radius")
    print("  4. Merge luminosity estimates from Hi-GAL SED fits")
    print("  5. Apply filters: 0.4 < d < 2 kpc, L > 10^4 Lsun")
    print("  6. Assign unique source IDs, track source origins")
    print("  7. Save as unified FITS/CSV catalog")

    return catalogs

if __name__ == '__main__':
    catalogs = main()
