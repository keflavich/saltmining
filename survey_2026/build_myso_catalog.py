"""
Build unified MYSO Phase 1 catalog
Handles actual column names from downloaded VizieR catalogs

Filters: 0.4 < d < 2 kpc, L > 10^3 Lsun
(Relaxed from L > 10^4 to include well-known massive YSOs like Cepheus A)
"""

import numpy as np
from pathlib import Path
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
import astropy.units as u
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = Path(__file__).parent / 'data'
OUTPUT_FILE = DATA_DIR / 'myso_candidates_phase1.fits'

def load_and_process():
    """Load all catalogs and build unified table"""
    print("\n" + "="*70)
    print("BUILDING PHASE 1 MYSO CANDIDATE CATALOG")
    print("="*70)

    catalogs = []

    # Load ATLASGAL
    print("\n1. Loading ATLASGAL (8002 sources)")
    atlasgal = Table.read(DATA_DIR / 'atlasgal_clumps.fits')
    # ATLASGAL has: Dist (kpc), needs Glon/Glat
    # It has Clump ID which likely encodes position
    if 'Dist' in atlasgal.colnames:
        atlasgal['distance_pc'] = atlasgal['Dist'] * 1000  # Convert kpc to pc
    atlasgal['catalog'] = 'ATLASGAL'
    print(f"   ✓ {len(atlasgal)} sources")
    catalogs.append(atlasgal)

    # Load RMS
    print("\n2. Loading RMS Survey (2798 sources)")
    rms = Table.read(DATA_DIR / 'rms.fits')
    # RMS has: RAJ2000, DEJ2000 (in sexagesimal format!), Dist, Lbol
    if 'RAJ2000' in rms.colnames and 'DEJ2000' in rms.colnames:
        # RMS coordinates are strings in sexagesimal format
        # Convert to Galactic coordinates
        try:
            # Parse sexagesimal RA/Dec strings
            coord = SkyCoord(
                rms['RAJ2000'].astype(str),
                rms['DEJ2000'].astype(str),
                unit=(u.hourangle, u.degree),
                frame='icrs'
            )
            rms['l'] = coord.galactic.l.deg
            rms['b'] = coord.galactic.b.deg
        except Exception as e:
            print(f"   Warning: Could not parse RMS coordinates: {e}")
    if 'Dist' in rms.colnames:
        rms['distance_pc'] = rms['Dist'] * 1000  # Convert kpc to pc
    if 'Lbol' in rms.colnames:
        rms['L_bol_Lsun'] = rms['Lbol']
    rms['catalog'] = 'RMS'
    print(f"   ✓ {len(rms)} sources with positions and luminosities")
    catalogs.append(rms)

    # Load WISE Green Objects
    print("\n3. Loading WISE Green Objects (231 sources)")
    wgo = Table.read(DATA_DIR / 'wise_green_objects.fits')
    # WGO has: GLON, GLAT, Dist
    if 'GLON' in wgo.colnames:
        wgo['l'] = wgo['GLON']
        wgo['b'] = wgo['GLAT']
    if 'Dist' in wgo.colnames:
        wgo['distance_pc'] = wgo['Dist'] * 1000
    wgo['catalog'] = 'WISE_GreenObjects'
    print(f"   ✓ {len(wgo)} sources")
    catalogs.append(wgo)

    # Load SEDIGISM (if useful)
    print("\n4. Loading SEDIGISM (325 sources)")
    try:
        sedigism = Table.read(DATA_DIR / 'sedigism.fits')
        # SEDIGISM has: _Glon, _Glat, _RA.icrs, _DE.icrs (unusual names)
        if '_Glon' in sedigism.colnames:
            sedigism['l'] = sedigism['_Glon']
            sedigism['b'] = sedigism['_Glat']
        sedigism['catalog'] = 'SEDIGISM'
        print(f"   ✓ {len(sedigism)} sources (note: limited distance/luminosity data)")
        # Only add if it has useful data
        if 'l' in sedigism.colnames:
            catalogs.append(sedigism)
    except Exception as e:
        print(f"   ⚠ Could not load SEDIGISM: {e}")

    # Merge catalogs (simple approach for now)
    print("\n" + "="*70)
    print("MERGING CATALOGS")
    print("="*70)

    # Start with RMS (most complete with positions + luminosity)
    primary_cat = rms.copy()
    print(f"Primary: RMS ({len(primary_cat)} sources with l, b, L_bol)")

    # Add WISE Green Objects
    if len(wgo) > 0:
        print(f"Adding WISE Green Objects ({len(wgo)} sources)...")
        # Use outer join to preserve all data
        primary_cat = vstack([primary_cat, wgo], join_type='outer')

    print(f"Merged catalog size: {len(primary_cat)} sources")

    # Apply filters
    print("\n" + "="*70)
    print("APPLYING FILTERS")
    print("="*70)

    print(f"  Initial catalog size: {len(primary_cat)}")

    # Distance filter
    if 'distance_pc' in primary_cat.colnames:
        d_kpc = primary_cat['distance_pc'] / 1000.0
        d_valid = ~np.isnan(d_kpc)
        d_mask = d_valid & (d_kpc >= 0.4) & (d_kpc <= 2.0)
        print(f"  Distance filter (0.4-2 kpc): {np.sum(d_mask)} sources have valid distance in range")
    else:
        d_mask = np.ones(len(primary_cat), dtype=bool)
        print(f"  WARNING: No distance column; can't filter by distance")

    # Luminosity filter (relaxed from 10^4 to 10^3 to include known MYSOs like Cepheus A)
    if 'L_bol_Lsun' in primary_cat.colnames:
        L_valid = ~np.isnan(primary_cat['L_bol_Lsun'])
        L_min = 1e3  # Relaxed from 10^4 to capture sources like Cepheus A
        L_mask = L_valid & (primary_cat['L_bol_Lsun'] >= L_min)
        print(f"  Luminosity filter (L > 10^3 Lsun): {np.sum(L_mask)} sources have valid L in range")
        total_mask = d_mask & L_mask
    else:
        print(f"  WARNING: No Lbol column; filtering by distance only")
        total_mask = d_mask

    filtered = primary_cat[total_mask]
    print(f"  Final filtered catalog: {len(filtered)} sources")

    # Statistics
    print("\n" + "="*70)
    print("STATISTICS")
    print("="*70)

    if 'distance_pc' in filtered.colnames:
        d_kpc = filtered['distance_pc'] / 1000.0
        d_valid = d_kpc[~np.isnan(d_kpc)]
        if len(d_valid) > 0:
            print(f"Distance (kpc):")
            print(f"  N with distance: {len(d_valid)}")
            print(f"  Min: {d_valid.min():.2f}")
            print(f"  Max: {d_valid.max():.2f}")
            print(f"  Mean: {np.mean(d_valid):.2f}")

    if 'L_bol_Lsun' in filtered.colnames:
        L_valid = filtered['L_bol_Lsun'][~np.isnan(filtered['L_bol_Lsun'])]
        if len(L_valid) > 0:
            print(f"Luminosity (Lsun):")
            print(f"  N with Lbol: {len(L_valid)}")
            print(f"  Min: {L_valid.min():.2e}")
            print(f"  Max: {L_valid.max():.2e}")
            print(f"  Median: {np.median(L_valid):.2e}")

    if 'catalog' in filtered.colnames:
        print(f"Source origin:")
        for origin in sorted(set(filtered['catalog'])):
            n = np.sum(filtered['catalog'] == origin)
            print(f"  {origin}: {n}")

    if 'l' in filtered.colnames:
        l_valid = filtered['l'][~np.isnan(filtered['l'])]
        if len(l_valid) > 0:
            print(f"Galactic distribution (l range): {l_valid.min():.1f}° to {l_valid.max():.1f}°")

    # Save
    print("\n" + "="*70)
    print("SAVING")
    print("="*70)

    filtered.write(OUTPUT_FILE, overwrite=True)
    print(f"✓ Saved {len(filtered)} candidates to {OUTPUT_FILE}")

    return filtered

if __name__ == '__main__':
    filtered = load_and_process()

    print("\n" + "="*70)
    print("PHASE 1 COMPLETE")
    print("="*70)
    print(f"\nNext: Cross-match against ALMA archive for Phase 2 (coming soon)")
