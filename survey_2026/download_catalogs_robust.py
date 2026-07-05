"""
Robust catalog download with multiple fallback strategies
Tries VizieR first, provides direct download links as fallback
"""

import urllib.request
import urllib.error
from pathlib import Path
from astroquery.vizier import Vizier
from astropy.table import Table
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = Path(__file__).parent / 'data'
DATA_DIR.mkdir(exist_ok=True)

# Catalog information with multiple ID attempts and fallback URLs
CATALOGS = {
    'atlasgal': {
        'description': 'ATLASGAL complete clump catalog (Urquhart 2018)',
        'vizier_ids': ['J/MNRAS/473/1059', '2018yCat..74731059U'],
        'downloaded': True,  # Already have this
    },
    'rms': {
        'description': 'Red MSX Source Survey (Lumsden+ 2013)',
        'vizier_ids': ['J/ApJS/208/11', '2013yCat..22080011L'],
        'fallback_url': 'http://rms.leeds.ac.uk/',
        'notes': 'Can download tables from RMS database at fallback_url'
    },
    'wise_green_objects': {
        'description': 'WISE Green Objects (Zhang et al. 2023)',
        'vizier_ids': ['J/ApJ/954/105', '2023yCat..22640024Z'],
        'notes': 'Try via direct ADS download if VizieR fails'
    },
    'higal': {
        'description': 'Hi-GAL SED fits (Berta+ 2013)',
        'vizier_ids': ['J/A+A/559/A19', '2013yCat..35510100B'],
    },
    'sedigism': {
        'description': 'SEDIGISM clouds (Schuller+ 2021)',
        'vizier_ids': ['J/MNRAS/500/3064', 'J/A+A/646/A47', '2021MNRAS.500.3064S'],
        'fallback_url': 'http://sedigism.mpifr-bonn.mpg.de/',
        'notes': 'Can download from SEDIGISM database browser'
    }
}

def download_vizier_catalog(catalog_name, vizier_ids, timeout=120):
    """Try to download catalog from VizieR with multiple ID attempts"""
    print(f"\n{'='*60}")
    print(f"Attempting: {catalog_name}")
    print(f"{'='*60}")

    v = Vizier(columns=['all'])
    v.ROW_LIMIT = -1
    v.timeout = timeout

    for vizier_id in vizier_ids:
        print(f"  Trying VizieR ID: {vizier_id}")
        try:
            result = v.query_constraints(catalog=vizier_id)
            if result and len(result) > 0:
                table = result[0]
                print(f"  ✓ SUCCESS: Downloaded {len(table)} sources")
                return table
            else:
                print(f"    No results for this ID")
        except Exception as e:
            print(f"    Failed: {type(e).__name__}")
            continue

    print(f"  ✗ All VizieR IDs failed for {catalog_name}")
    return None

def save_catalog(table, catalog_name):
    """Save table to FITS file"""
    if table is None:
        return False

    filename = f"{catalog_name.lower()}.fits"
    output_path = DATA_DIR / filename

    try:
        table.write(output_path, overwrite=True)
        print(f"  ✓ Saved {len(table)} sources to {output_path}")
        return True
    except Exception as e:
        print(f"  ✗ Error saving: {e}")
        return False

def print_fallback_info(catalog_name, info):
    """Print fallback download information"""
    if 'fallback_url' in info:
        print(f"\n  FALLBACK: {info['fallback_url']}")
    if 'notes' in info:
        print(f"  Note: {info['notes']}")

def main():
    print("\n" + "="*70)
    print("ROBUST MYSO CATALOG DOWNLOAD")
    print("="*70)

    results = {}

    for catalog_name, catalog_info in CATALOGS.items():
        # Skip ATLASGAL (already downloaded)
        if catalog_info.get('downloaded'):
            print(f"\n✓ {catalog_name}: Already downloaded, skipping")
            continue

        print(f"\n{'='*70}")
        print(f"{catalog_info['description']}")
        print(f"{'='*70}")

        # Try VizieR
        table = download_vizier_catalog(
            catalog_name,
            catalog_info['vizier_ids'],
            timeout=120
        )

        if table is not None:
            if save_catalog(table, catalog_name):
                results[catalog_name] = 'SUCCESS'
            else:
                results[catalog_name] = 'DOWNLOAD_OK_SAVE_FAILED'
                print_fallback_info(catalog_name, catalog_info)
        else:
            results[catalog_name] = 'VIZIER_FAILED'
            print_fallback_info(catalog_name, catalog_info)

    # Summary
    print("\n" + "="*70)
    print("DOWNLOAD SUMMARY")
    print("="*70)

    for catalog_name, status in sorted(results.items()):
        if status == 'SUCCESS':
            print(f"  ✓ {catalog_name:25} {status}")
        else:
            print(f"  ⚠ {catalog_name:25} {status}")

    print("\n" + "="*70)
    print("NEXT STEPS")
    print("="*70)
    print("""
For any failed catalogs, download manually from:
  - RMS: http://rms.leeds.ac.uk/ → download main source table
  - WISE Green Objects: https://ui.adsabs.harvard.edu/abs/2023yCat..22640024Z/abstract
  - Hi-GAL: https://ui.adsabs.harvard.edu/abs/2013yCat..35510100B/abstract
  - SEDIGISM: http://sedigism.mpifr-bonn.mpg.de/ (or arXiv for paper)

Place downloaded FITS files in: {DATA_DIR}

Then run: python load_and_crossmatch.py
""")

if __name__ == '__main__':
    main()
