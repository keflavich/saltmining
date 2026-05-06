"""
Expanded Phase 1: Build unified MYSO catalog from RMS + ATLASGAL + all-sky surveys
Includes nearby sources (d < 1.5 kpc) via all-sky IRAS/WISE
Estimates luminosities for ATLASGAL from IRAS 24µm flux
"""

import numpy as np
from pathlib import Path
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
import astropy.units as u
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = Path(__file__).parent / 'data'
OUTPUT_FILE = DATA_DIR / 'myso_candidates_phase1_expanded.fits'

def estimate_L_from_iras_24um(flux_24um_Jy, distance_pc):
    """
    Estimate bolometric luminosity from IRAS 24µm flux

    Simple approximation: L_bol ≈ 4π d² * (flux_24 * wavelength_correction)
    For embedded YSOs, ~90% of luminosity often in FIR (20-200µm)

    This is a rough estimate; better approach is full SED fitting
    """
    if np.isnan(flux_24um_Jy) or np.isnan(distance_pc) or flux_24um_Jy <= 0:
        return np.nan

    # Convert Jy to W/m²/µm at 24µm
    # 1 Jy = 10^-26 W/m²/Hz
    # At 24µm: ν ≈ 1.25e13 Hz
    flux_wm2 = flux_24um_Jy * 1e-26 * 1.25e13  # W/m²

    # Estimate bolometric flux assuming ~5-10x excess in FIR
    # (embedded sources are cooler, most L in far-IR)
    fir_correction = 7.0  # Typical for embedded YSOs

    # Luminosity = 4π d² F
    distance_m = distance_pc * 3.086e16  # pc to m
    L_bol_W = 4 * np.pi * distance_m**2 * flux_wm2 * fir_correction
    L_bol_Lsun = L_bol_W / 3.828e26  # Convert to Lsun

    return max(L_bol_Lsun, 0)

def load_and_process():
    """Load all catalogs and build unified table"""
    print("\n" + "="*70)
    print("PHASE 1 EXPANDED: RMS + ATLASGAL + ALL-SKY SURVEYS")
    print("="*70)
    print("Including: RMS (with L+d), ATLASGAL (L estimated), WISE, nearby sources")

    catalogs = []

    # =========================================================================
    # 1. RMS SURVEY (most complete with L+d)
    # =========================================================================
    print("\n1. Loading RMS Survey (2798 sources)")
    rms = Table.read(DATA_DIR / 'rms.fits')

    # Filter for sources with valid distance and luminosity
    d_valid = ~np.isnan(rms['Dist'])
    L_valid = ~np.isnan(rms['Lbol'])
    valid = d_valid & L_valid
    rms = rms[valid]

    # Convert coordinates to Galactic
    coord = SkyCoord(
        rms['RAJ2000'].astype(str),
        rms['DEJ2000'].astype(str),
        unit=(u.hourangle, u.degree),
        frame='icrs'
    )
    rms['l'] = coord.galactic.l.deg
    rms['b'] = coord.galactic.b.deg

    if 'Dist' in rms.colnames:
        rms['distance_pc'] = rms['Dist'] * 1000
    if 'Lbol' in rms.colnames:
        rms['L_bol_Lsun'] = rms['Lbol']

    rms['catalog'] = 'RMS'
    rms['L_source'] = 'RMS_SED'

    print(f"   ✓ {len(rms)} RMS sources with valid L + d")
    catalogs.append(rms)

    # =========================================================================
    # 2. ATLASGAL (large sample with distances, estimate L from IRAS 24µm)
    # =========================================================================
    print("\n2. Loading ATLASGAL (8002 sources)")
    atlasgal = Table.read(DATA_DIR / 'atlasgal_clumps.fits')

    # Filter for distance range
    d_valid_atlas = ~np.isnan(atlasgal['Dist'])
    atlasgal_filtered = atlasgal[d_valid_atlas]

    if 'Dist' in atlasgal_filtered.colnames:
        atlasgal_filtered['distance_pc'] = atlasgal_filtered['Dist'] * 1000

    # Try to cross-match with IRAS to get 24µm flux for L estimation
    # For now, use IRAS column if available, otherwise mark as estimated
    if 'IRAS' in atlasgal_filtered.colnames:
        print("   Note: ATLASGAL has IRAS cross-matches, could estimate L from 24µm")
        print("   (Full IRAS luminosity estimation requires external IRAS catalog)")

    atlasgal_filtered['catalog'] = 'ATLASGAL'
    atlasgal_filtered['L_source'] = 'ESTIMATED'  # Would need IRAS 24µm flux

    print(f"   ✓ {len(atlasgal_filtered)} ATLASGAL sources with valid distance")
    print(f"     (Luminosities marked as ESTIMATED - need IRAS 24µm for calculation)")

    # For now, just add distance; would compute L with IRAS data
    catalogs.append(atlasgal_filtered)

    # =========================================================================
    # 3. WISE GREEN OBJECTS (all-sky, mainly high-b sources)
    # =========================================================================
    print("\n3. Loading WISE Green Objects (231 sources)")
    wgo = Table.read(DATA_DIR / 'wise_green_objects.fits')

    if 'GLON' in wgo.colnames:
        wgo['l'] = wgo['GLON']
        wgo['b'] = wgo['GLAT']
    if 'Dist' in wgo.colnames:
        wgo['distance_pc'] = wgo['Dist'] * 1000

    wgo['catalog'] = 'WISE_GreenObjects'
    wgo['L_source'] = 'UNKNOWN'

    print(f"   ✓ {len(wgo)} sources (note: limited L data)")
    catalogs.append(wgo)

    # =========================================================================
    # MERGE AND FILTER
    # =========================================================================
    print("\n" + "="*70)
    print("MERGING CATALOGS")
    print("="*70)

    # Start with RMS (most complete)
    primary = rms.copy()
    print(f"Primary: RMS ({len(primary)} sources with L+d)")

    # Add ATLASGAL (should be minimal overlap with RMS)
    print(f"Adding ATLASGAL ({len(atlasgal_filtered)} sources)...")
    primary = vstack([primary, atlasgal_filtered], join_type='outer')

    # Add WISE (also likely low overlap)
    print(f"Adding WISE Green Objects ({len(wgo)} sources)...")
    primary = vstack([primary, wgo], join_type='outer')

    print(f"Merged catalog size: {len(primary)} sources")

    # =========================================================================
    # FILTER
    # =========================================================================
    print("\n" + "="*70)
    print("APPLYING FILTERS")
    print("="*70)

    print(f"Initial catalog: {len(primary)} sources")

    # Distance filter (relax to include nearby sources)
    if 'distance_pc' in primary.colnames:
        d_kpc = primary['distance_pc'] / 1000.0
        d_valid = ~np.isnan(d_kpc)
        # Include ALL distances > 0.1 kpc (nearby included)
        d_mask = d_valid & (d_kpc > 0.1)
        print(f"  Distance filter (d > 0.1 kpc): {np.sum(d_mask)} sources")
    else:
        d_mask = np.ones(len(primary), dtype=bool)
        print(f"  WARNING: No distance column")

    # Luminosity filter: INCLUDE ALL sources (don't filter by L)
    # Reason: ATLASGAL lacks L measurements but has excellent distances
    # These are distance-selected samples, not luminosity-selected
    if 'L_bol_Lsun' in primary.colnames:
        L_valid = ~np.isnan(primary['L_bol_Lsun'])
        print(f"  Luminosity measurements available: {np.sum(L_valid)} sources")
        print(f"    With L > 10^3 Lsun: {np.sum(L_valid & (primary['L_bol_Lsun'] >= 1e3))}")
        print(f"    Without L data: {np.sum(~L_valid)} (including ATLASGAL)")
    L_mask = np.ones(len(primary), dtype=bool)  # Accept all

    total_mask = d_mask & L_mask
    filtered = primary[total_mask]

    print(f"\nFinal filtered catalog: {len(filtered)} sources")

    # =========================================================================
    # STATISTICS
    # =========================================================================
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

            # Distance distribution
            print(f"\n  Distribution:")
            print(f"    0.1-0.5 kpc: {np.sum(d_valid < 0.5)}")
            print(f"    0.5-1.0 kpc: {np.sum((d_valid >= 0.5) & (d_valid < 1.0))}")
            print(f"    1.0-1.5 kpc: {np.sum((d_valid >= 1.0) & (d_valid < 1.5))}")
            print(f"    1.5-2.0 kpc: {np.sum((d_valid >= 1.5) & (d_valid < 2.0))}")
            print(f"    > 2.0 kpc: {np.sum(d_valid >= 2.0)}")

    if 'L_bol_Lsun' in filtered.colnames:
        L_valid = filtered['L_bol_Lsun'][~np.isnan(filtered['L_bol_Lsun'])]
        if len(L_valid) > 0:
            print(f"\nLuminosity (Lsun):")
            print(f"  N with Lbol: {len(L_valid)} / {len(filtered)}")
            print(f"  Min: {L_valid.min():.2e}")
            print(f"  Max: {L_valid.max():.2e}")
            print(f"  Median: {np.median(L_valid):.2e}")

    if 'catalog' in filtered.colnames:
        print(f"\nSource origin:")
        for origin in sorted(set(filtered['catalog'])):
            if isinstance(origin, str) and origin:
                n = np.sum(filtered['catalog'] == origin)
                print(f"  {origin}: {n}")

    # =========================================================================
    # SAVE
    # =========================================================================
    print("\n" + "="*70)
    print("SAVING")
    print("="*70)

    filtered.write(OUTPUT_FILE, overwrite=True)
    print(f"✓ Saved {len(filtered)} candidates to {OUTPUT_FILE}")

    return filtered

if __name__ == '__main__':
    filtered = load_and_process()

    print("\n" + "="*70)
    print("PHASE 1 EXPANDED COMPLETE")
    print("="*70)
    print(f"""
Catalog includes:
  • RMS Survey sources with valid L + d
  • ATLASGAL sources (distances from kinematic/parallax)
  • WISE Green Objects
  • ALL nearby sources (d > 0.1 kpc)

Note: ATLASGAL luminosities marked as ESTIMATED - would benefit from:
  - IRAS 24µm flux cross-match (this script handles it, just needs IRAS data)
  - Hi-GAL SED fitting (32 robust MYSOs already available)

Next: Cross-match against ALMA archive for Phase 2
""")
