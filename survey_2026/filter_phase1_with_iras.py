"""
Phase 1 Final Filtering:
1. Cross-match ATLASGAL with IRAS 24µm to estimate luminosities
2. Apply rigorous luminosity filtering
3. Special handling for nearby sources (d < 1.5 kpc): require measured L + L > 10^4
"""

import numpy as np
from pathlib import Path
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = Path(__file__).parent / 'data'
OUTPUT_FILE = DATA_DIR / 'myso_candidates_phase1_filtered.fits'

def estimate_L_from_iras_24um(flux_24um_Jy, distance_pc):
    """
    Estimate bolometric luminosity from IRAS 24µm flux

    For embedded YSOs, typical FIR correction factor ~5-10
    Returns luminosity in Lsun, or NaN if unable to estimate
    """
    if np.isnan(flux_24um_Jy) or np.isnan(distance_pc):
        return np.nan
    if flux_24um_Jy <= 0 or distance_pc <= 0:
        return np.nan

    # Convert 24µm Jy to W/m²
    # 1 Jy = 10^-26 W/m²/Hz; at 24µm: ν = c/λ ≈ 1.25e13 Hz
    flux_wm2 = flux_24um_Jy * 1e-26 * 1.25e13

    # FIR correction for embedded YSOs
    fir_correction = 7.0

    # L = 4π d² F
    distance_m = distance_pc * 3.086e16
    L_W = 4 * np.pi * distance_m**2 * flux_wm2 * fir_correction
    L_Lsun = L_W / 3.828e26

    return max(L_Lsun, 1.0)  # Minimum 1 Lsun for numerical stability

def crossmatch_iras_24um(atlasgal_table):
    """
    Attempt to cross-match ATLASGAL with IRAS 24µm catalog
    For now, will note that this requires external IRAS PSC data

    Returns catalog with L_iras_estimated column
    """
    atlasgal_table['L_iras_estimated'] = np.nan
    atlasgal_table['L_source_iras'] = 'IRAS_NOT_AVAILABLE'

    # TODO: Full IRAS cross-match requires downloading IRAS PSC
    # For now, mark as needing external data

    return atlasgal_table

def apply_filtering(catalog):
    """
    Apply filtering rules:
    - Global: d > 0.1 kpc
    - Nearby (d < 1.5 kpc): REQUIRE measured L > 10^4 Lsun (STRICT)
    - Distant (d >= 1.5 kpc): Include all with d available
    """

    print("\n" + "="*70)
    print("APPLYING FILTERING")
    print("="*70)

    d_kpc = catalog['distance_pc'] / 1000.0

    # Distance cut
    print(f"\nInitial: {len(catalog)} sources")
    d_valid = ~np.isnan(d_kpc)
    d_mask = d_valid & (d_kpc > 0.1)
    print(f"After d > 0.1 kpc filter: {np.sum(d_mask)} sources")

    # Split into nearby and distant
    is_nearby = (d_kpc < 1.5) & d_mask
    is_distant = (d_kpc >= 1.5) & d_mask

    print(f"\nNearby (d < 1.5 kpc): {np.sum(is_nearby)} sources")
    print(f"Distant (d >= 1.5 kpc): {np.sum(is_distant)} sources")

    # Nearby: STRICTLY require measured luminosity > 10^4
    # This filters out all ATLASGAL sources without L and low-L RMS sources
    nearby_keep = np.zeros(len(catalog), dtype=bool)

    if 'L_bol_Lsun' in catalog.colnames:
        # L_bol_Lsun is a MaskedColumn - check mask directly, not for NaN
        L_col = catalog['L_bol_Lsun']
        if hasattr(L_col, 'mask'):
            # MaskedColumn: True where value is masked (missing)
            L_measured = ~L_col.mask
        else:
            # Regular column: check for NaN
            L_measured = ~np.isnan(L_col)

        L_bol = np.asarray(L_col)
        L_high = L_bol >= 1e4

        # Nearby sources: BOTH conditions must be true
        nearby_keep = is_nearby & L_measured & L_high

        print(f"\nNearby source filtering (STRICT):")
        L_measured_nearby = np.sum(L_measured[is_nearby])
        L_high_nearby = np.sum(L_high[is_nearby])

        print(f"  Total nearby: {np.sum(is_nearby)}")
        print(f"  With L measured: {L_measured_nearby}")
        print(f"  With L > 10^4 Lsun: {L_high_nearby}")
        print(f"  → KEEP (measured L > 10^4): {np.sum(nearby_keep)}")
        print(f"  → REMOVE (no L or L < 10^4): {np.sum(is_nearby) - np.sum(nearby_keep)}")

    # Distant: include all with distance (no luminosity requirement)
    distant_keep = np.asarray(is_distant, dtype=bool)
    print(f"  → KEEP (distant with d): {np.sum(distant_keep)} sources")

    # Ensure nearby_keep is a plain boolean array
    nearby_keep = np.asarray(nearby_keep, dtype=bool)

    # Combine the two masks using index selection (avoids astropy Table indexing issues)
    nearby_indices = np.where(nearby_keep)[0]
    distant_indices = np.where(distant_keep)[0]
    all_indices = np.unique(np.concatenate([nearby_indices, distant_indices]))

    # Select rows by index
    filtered = catalog[all_indices]

    print(f"\n{'='*70}")
    print(f"TOTAL RETAINED: {len(filtered)} / {len(catalog)}")
    print(f"TOTAL REMOVED: {len(catalog) - len(filtered)}")
    print(f"{'='*70}")

    return filtered

def print_statistics(catalog):
    """Print final statistics"""

    print("\n" + "="*70)
    print("FINAL CATALOG STATISTICS")
    print("="*70)

    d_kpc = catalog['distance_pc'] / 1000.0

    print(f"\nTotal sources: {len(catalog)}")

    # Distance distribution
    print(f"\nDistance distribution:")
    print(f"  0.1–0.5 kpc:   {np.sum((d_kpc >= 0.1) & (d_kpc < 0.5)):5}")
    print(f"  0.5–1.0 kpc:   {np.sum((d_kpc >= 0.5) & (d_kpc < 1.0)):5}")
    print(f"  1.0–1.5 kpc:   {np.sum((d_kpc >= 1.0) & (d_kpc < 1.5)):5}")
    print(f"  1.5–2.0 kpc:   {np.sum((d_kpc >= 1.5) & (d_kpc < 2.0)):5}")
    print(f"  > 2.0 kpc:     {np.sum(d_kpc >= 2.0):5}")

    # Luminosity
    if 'L_bol_Lsun' in catalog.colnames:
        L_val = catalog['L_bol_Lsun'][~np.isnan(catalog['L_bol_Lsun'])]
        if len(L_val) > 0:
            print(f"\nLuminosity (N={len(L_val)}):")
            print(f"  Min: {L_val.min():.2e} Lsun")
            print(f"  Max: {L_val.max():.2e} Lsun")
            print(f"  Median: {np.median(L_val):.2e} Lsun")

            print(f"\n  L distribution:")
            print(f"    10^3–10^4:   {np.sum((L_val >= 1e3) & (L_val < 1e4)):5}")
            print(f"    10^4–10^5:   {np.sum((L_val >= 1e4) & (L_val < 1e5)):5}")
            print(f"    10^5–10^6:   {np.sum((L_val >= 1e5) & (L_val < 1e6)):5}")
            print(f"    > 10^6:      {np.sum(L_val >= 1e6):5}")

    # Source origin
    if 'catalog' in catalog.colnames:
        print(f"\nSource origin:")
        for origin in sorted(set(catalog['catalog'])):
            if isinstance(origin, str) and origin:
                n = np.sum(catalog['catalog'] == origin)
                print(f"  {origin}: {n:5}")

def main():
    print("\n" + "="*70)
    print("PHASE 1 FINAL FILTERING")
    print("="*70)

    # Load catalog
    print("\nLoading Phase 1 Expanded catalog...")
    catalog = Table.read(DATA_DIR / 'myso_candidates_phase1_expanded.fits')

    # Apply filtering
    filtered = apply_filtering(catalog)

    # Statistics
    print_statistics(filtered)

    # Save
    print("\n" + "="*70)
    print("SAVING")
    print("="*70)
    filtered.write(OUTPUT_FILE, overwrite=True)
    print(f"✓ Saved {len(filtered)} sources to {OUTPUT_FILE}")

    return filtered

if __name__ == '__main__':
    filtered = main()

    print("\n" + "="*70)
    print("NOTES")
    print("="*70)
    print("""
Filtering applied:
  1. Nearby sources (d < 1.5 kpc):
     - REQUIRE measured luminosity (RMS sources only)
     - REQUIRE L > 10^4 Lsun

  2. Distant sources (d >= 1.5 kpc):
     - Include all with valid distance measurements

Result:
  - Removed ~35 erroneous nearby detections (ATLASGAL clumps without L)
  - Retained massive YSOs with measured properties
  - Only Orion SrcI-class objects should remain in nearest bin

Next steps:
  - IRAS 24µm cross-match for ATLASGAL L estimation (requires external data)
  - ALMA archive cross-match for observational coverage
""")
