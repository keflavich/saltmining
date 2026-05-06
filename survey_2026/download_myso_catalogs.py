"""
Download MYSO catalogs directly from VizieR using astroquery
Phase 1: Build candidate MYSO table with L > 10^4 Lsun, 0.4 < d < 2 kpc
"""

import numpy as np
from pathlib import Path
import astropy.units as u
from astropy.table import Table, vstack
from astroquery.vizier import Vizier
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = Path(__file__).parent / 'data'
DATA_DIR.mkdir(exist_ok=True)

def download_atlasgal():
    """
    Download ATLASGAL complete clump catalog with distances (Urquhart 2018)
    VizieR catalog: J/MNRAS/473/1059
    """
    print("\n" + "="*60)
    print("DOWNLOADING: ATLASGAL (Urquhart 2018)")
    print("="*60)
    try:
        v = Vizier(columns=['all'])
        v.ROW_LIMIT = -1

        # Query ATLASGAL main table
        result = v.query_constraints(catalog='J/MNRAS/473/1059')

        if not result:
            print("ERROR: VizieR query returned no results")
            return None

        atlasgal = result[0]
        print(f"✓ Downloaded {len(atlasgal)} ATLASGAL sources")
        print(f"  Columns: {atlasgal.colnames[:15]}...")

        # Save to FITS
        output_path = DATA_DIR / 'atlasgal_clumps.fits'
        atlasgal.write(output_path, overwrite=True)
        print(f"✓ Saved to {output_path}")
        return atlasgal

    except Exception as e:
        print(f"✗ Error downloading ATLASGAL: {e}")
        return None

def download_rms_survey():
    """
    Download Red MSX Source (RMS) Survey catalog (Lumsden et al. 2013)
    VizieR catalog: J/ApJS/208/11
    ~1200 massive YSOs and UCHII regions
    """
    print("\n" + "="*60)
    print("DOWNLOADING: RMS Survey (Lumsden et al. 2013)")
    print("="*60)
    try:
        v = Vizier(columns=['all'])
        v.ROW_LIMIT = -1

        # Query RMS main catalog (EVLA follow-up observations)
        result = v.query_constraints(catalog='J/ApJS/208/11')

        if not result:
            print("ERROR: VizieR query returned no results")
            return None

        rms = result[0]
        print(f"✓ Downloaded {len(rms)} RMS sources")
        print(f"  Columns: {rms.colnames[:15]}...")

        # Save to FITS
        output_path = DATA_DIR / 'rms_survey.fits'
        rms.write(output_path, overwrite=True)
        print(f"✓ Saved to {output_path}")
        return rms

    except Exception as e:
        print(f"✗ Error downloading RMS: {e}")
        return None

def download_wise_green_objects():
    """
    Download WISE Green Objects catalog (Zhang et al. 2023)
    VizieR catalog: J/ApJ/954/105
    """
    print("\n" + "="*60)
    print("DOWNLOADING: WISE Green Objects (Zhang et al. 2023)")
    print("="*60)
    try:
        v = Vizier(columns=['all'])
        v.ROW_LIMIT = -1

        result = v.query_constraints(catalog='J/ApJ/954/105')

        if not result:
            print("ERROR: VizieR query returned no results")
            return None

        wgo = result[0]
        print(f"✓ Downloaded {len(wgo)} WISE Green Objects")
        print(f"  Columns: {wgo.colnames[:15]}...")

        # Save to FITS
        output_path = DATA_DIR / 'wise_green_objects.fits'
        wgo.write(output_path, overwrite=True)
        print(f"✓ Saved to {output_path}")
        return wgo

    except Exception as e:
        print(f"✗ Error downloading WISE Green Objects: {e}")
        return None

def download_higal_sources():
    """
    Download Hi-GAL compact source catalog with SED-fitted luminosities
    Multiple tables available; download Berta et al. 2013 panchromatic SED
    VizieR catalog: J/A+A/559/A19
    """
    print("\n" + "="*60)
    print("DOWNLOADING: Hi-GAL SED Fits (Berta et al. 2013)")
    print("="*60)
    try:
        v = Vizier(columns=['all'])
        v.ROW_LIMIT = -1

        result = v.query_constraints(catalog='J/A+A/559/A19')

        if not result:
            print("ERROR: VizieR query returned no results")
            return None

        higal = result[0]
        print(f"✓ Downloaded {len(higal)} Hi-GAL sources")
        print(f"  Columns: {higal.colnames[:15]}...")

        # Save to FITS
        output_path = DATA_DIR / 'higal_sed_fits.fits'
        higal.write(output_path, overwrite=True)
        print(f"✓ Saved to {output_path}")
        return higal

    except Exception as e:
        print(f"✗ Error downloading Hi-GAL: {e}")
        return None

def download_sedigism_via_vizier():
    """
    Download SEDIGISM cloud catalog via VizieR
    VizieR catalog: J/A+A/646/A47 (Duarte-Cabral et al. 2021) or similar
    Note: SEDIGISM full catalog may be very large; check what's available
    """
    print("\n" + "="*60)
    print("DOWNLOADING: SEDIGISM clouds")
    print("="*60)
    try:
        v = Vizier(columns=['all'])
        v.ROW_LIMIT = -1

        # Try Duarte-Cabral et al. 2021 filaments catalog
        result = v.query_constraints(catalog='J/A+A/646/A47')

        if result:
            sedigism = result[0]
            print(f"✓ Downloaded {len(sedigism)} SEDIGISM sources")
            print(f"  Columns: {sedigism.colnames[:15]}...")

            output_path = DATA_DIR / 'sedigism_clouds.fits'
            sedigism.write(output_path, overwrite=True)
            print(f"✓ Saved to {output_path}")
            return sedigism
        else:
            print("⚠ SEDIGISM: No results from VizieR query")
            print("  NOTE: SEDIGISM main cloud catalog may need direct download from:")
            print("  http://sedigism.mpifr-bonn.mpg.de/")
            return None

    except Exception as e:
        print(f"⚠ Warning querying SEDIGISM: {e}")
        print("  Try downloading directly from SEDIGISM database")
        return None

def main():
    print("\n" + "="*70)
    print("DOWNLOADING MYSO CATALOGS - PHASE 1")
    print("="*70)

    catalogs = {}

    # Download in order
    catalogs['ATLASGAL'] = download_atlasgal()
    catalogs['RMS'] = download_rms_survey()
    catalogs['WISE_GreenObjects'] = download_wise_green_objects()
    catalogs['HiGAL'] = download_higal_sources()
    catalogs['SEDIGISM'] = download_sedigism_via_vizier()

    # Summary
    print("\n" + "="*70)
    print("DOWNLOAD SUMMARY")
    print("="*70)

    for name, cat in catalogs.items():
        if cat is not None:
            print(f"✓ {name:20} {len(cat):6} sources")
        else:
            print(f"✗ {name:20} FAILED or not available")

    print("\n" + "="*70)
    print("NEXT STEP")
    print("="*70)
    print("\nRun the cross-matching pipeline:")
    print("  python load_and_crossmatch.py")
    print("\nThis will:")
    print("  - Load all downloaded catalogs")
    print("  - Cross-match by Galactic position")
    print("  - Use Edenhofer dust map for distance verification")
    print("  - Apply filters: 0.4 < d < 2 kpc, L > 10^4 Lsun")
    print("  - Output: myso_candidates_phase1.fits")

if __name__ == '__main__':
    main()
