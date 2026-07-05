"""
Load and cross-match MYSO catalogs from downloaded FITS/ASCII files
Phase 1: Build unified candidate MYSO table
Includes ATLASGAL, RMS Survey, SEDIGISM, WISE Green Objects, Hi-GAL
Uses Edenhofer dust map for distance verification
"""

import numpy as np
import os
from pathlib import Path
import astropy.units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.table import Table, vstack, unique, Column, hstack
from astropy.io import fits
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = Path(__file__).parent / 'data'
DATA_DIR.mkdir(exist_ok=True)
DUSTMAP_PATH = Path('/orange/adamginsburg/galactic_plane_surveys/dustmaps/edenhofer_mean_and_std_lbd.fits')

class MYSOCatalogMerger:
    def __init__(self, data_dir=None):
        """Initialize catalog merger"""
        self.data_dir = data_dir or DATA_DIR
        self.catalogs = {}
        self.merged = None

    def load_atlasgal(self, filename='atlasgal_clumps.fits'):
        """Load ATLASGAL clump catalog with distances"""
        filepath = self.data_dir / filename
        if not filepath.exists():
            print(f"⚠ {filepath} not found. Download from:")
            print("  https://ui.adsabs.harvard.edu/abs/2018yCat..74731059U/abstract")
            return None

        print(f"Loading ATLASGAL from {filepath}")
        atlasgal = Table.read(filepath)

        # Standardize column names
        self._standardize_columns(atlasgal, {
            'GLON': 'l', 'Glon': 'l', 'glon': 'l',
            'GLAT': 'b', 'Glat': 'b', 'glat': 'b',
            'Dist': 'distance_pc', 'dist': 'distance_pc',
            'Lbol': 'L_bol_Lsun', 'Lbol_bol': 'L_bol_Lsun'
        })

        atlasgal['catalog_origin'] = 'ATLASGAL'
        self.catalogs['ATLASGAL'] = atlasgal
        print(f"  Loaded {len(atlasgal)} sources")
        return atlasgal

    def load_sedigism(self, filename='sedigism.fits'):
        """Load SEDIGISM cloud/clump catalog"""
        filepath = self.data_dir / filename
        if not filepath.exists():
            print(f"⚠ {filepath} not found. Download from:")
            print("  http://sedigism.mpifr-bonn.mpg.de/ or")
            print("  https://ui.adsabs.harvard.edu/abs/2021MNRAS.500.3064S/abstract")
            return None

        print(f"Loading SEDIGISM from {filepath}")
        sedigism = Table.read(filepath)

        self._standardize_columns(sedigism, {
            'GLON': 'l', 'lon': 'l',
            'GLAT': 'b', 'lat': 'b',
            'Dist': 'distance_pc', 'distance': 'distance_pc',
            'Lbol': 'L_bol_Lsun', 'Lum': 'L_bol_Lsun'
        })

        sedigism['catalog_origin'] = 'SEDIGISM'
        self.catalogs['SEDIGISM'] = sedigism
        print(f"  Loaded {len(sedigism)} sources")
        return sedigism

    def load_wise_green_objects(self, filename='wise_green_objects.fits'):
        """Load WISE Green Objects catalog - 231 robust MYSOs"""
        filepath = self.data_dir / filename
        if not filepath.exists():
            print(f"⚠ {filepath} not found. Download from:")
            print("  https://ui.adsabs.harvard.edu/abs/2023yCat..22640024Z/abstract")
            return None

        print(f"Loading WISE Green Objects from {filepath}")
        wgo = Table.read(filepath)

        # WGO has RA/Dec; convert to Galactic
        if 'RAJ2000' in wgo.colnames and 'DEJ2000' in wgo.colnames:
            coord = SkyCoord(ra=wgo['RAJ2000'], dec=wgo['DEJ2000'], frame='icrs')
            wgo['l'] = coord.galactic.l.deg
            wgo['b'] = coord.galactic.b.deg
        elif 'RA_2000' in wgo.colnames:
            coord = SkyCoord(ra=wgo['RA_2000'], dec=wgo['DE_2000'], frame='icrs')
            wgo['l'] = coord.galactic.l.deg
            wgo['b'] = coord.galactic.b.deg

        # Note: WGO needs cross-match with Hi-GAL for luminosity
        wgo['catalog_origin'] = 'WISE_GreenObjects'
        self.catalogs['WGO'] = wgo
        print(f"  Loaded {len(wgo)} sources (luminosity from Hi-GAL cross-match)")
        return wgo

    def load_higal(self, filename='higal.fits'):
        """Load Hi-GAL sources with SED-fitted luminosities"""
        filepath = self.data_dir / filename
        if not filepath.exists():
            print(f"⚠ {filepath} not found. Download from:")
            print("  https://ui.adsabs.harvard.edu/abs/2013yCat..35510100B/abstract")
            return None

        print(f"Loading Hi-GAL from {filepath}")
        higal = Table.read(filepath)

        # Convert to Galactic if needed
        if 'RAJ2000' in higal.colnames:
            coord = SkyCoord(ra=higal['RAJ2000'], dec=higal['DEJ2000'], frame='icrs')
            higal['l'] = coord.galactic.l.deg
            higal['b'] = coord.galactic.b.deg

        self._standardize_columns(higal, {
            'Lbol': 'L_bol_Lsun', 'Lbol_fit': 'L_bol_Lsun'
        })

        higal['catalog_origin'] = 'HiGAL'
        self.catalogs['HiGAL'] = higal
        print(f"  Loaded {len(higal)} sources")
        return higal

    def load_rms_survey(self, filename='rms.fits'):
        """Load Red MSX Source (RMS) Survey catalog (Lumsden et al. 2013)
        ~2800 massive YSOs and UCHII regions"""
        filepath = self.data_dir / filename
        if not filepath.exists():
            print(f"⚠ {filepath} not found. Download from:")
            print("  https://ui.adsabs.harvard.edu/abs/2013yCat..22080011L/abstract")
            return None

        print(f"Loading RMS Survey from {filepath}")
        rms = Table.read(filepath)

        # RMS typically has RA/Dec or Glon/Glat; convert as needed
        if 'RAJ2000' in rms.colnames and 'DEJ2000' in rms.colnames:
            coord = SkyCoord(ra=rms['RAJ2000'], dec=rms['DEJ2000'], frame='icrs')
            rms['l'] = coord.galactic.l.deg
            rms['b'] = coord.galactic.b.deg
        elif 'GLON' not in rms.colnames and 'Glon' in rms.colnames:
            rms['l'] = rms['Glon']
            rms['b'] = rms['Glat']

        self._standardize_columns(rms, {
            'GLON': 'l', 'Glon': 'l',
            'GLAT': 'b', 'Glat': 'b',
            'Dist': 'distance_pc', 'dist': 'distance_pc'
        })

        rms['catalog_origin'] = 'RMS'
        self.catalogs['RMS'] = rms
        print(f"  Loaded {len(rms)} sources")
        return rms

    @staticmethod
    def _standardize_columns(table, mapping):
        """Rename columns in table using mapping dict"""
        for old_names, new_name in mapping.items():
            if isinstance(old_names, str):
                old_names = [old_names]
            for old_name in old_names:
                if old_name in table.colnames:
                    table.rename_column(old_name, new_name)
                    break

    def load_edenhofer_dust_map(self):
        """Load Edenhofer & Zucker 3D dust extinction map
        Provides distances to clouds within 2 kpc"""
        print("\n=== Loading Edenhofer dust map ===")
        if not DUSTMAP_PATH.exists():
            print(f"⚠ Dust map not found at {DUSTMAP_PATH}")
            return None

        try:
            # Load the FITS file
            hdul = fits.open(DUSTMAP_PATH)
            dust_data = hdul[0].data
            dust_header = hdul[0].header
            hdul.close()

            print(f"✓ Loaded dust map from {DUSTMAP_PATH}")
            print(f"  Shape: {dust_data.shape}")
            print(f"  (This is a 3D map: l, b, distance)")
            return dust_data, dust_header
        except Exception as e:
            print(f"⚠ Error loading dust map: {e}")
            return None

    def crossmatch_by_position(self, match_radius=5*u.arcsec):
        """
        Cross-match catalogs by Galactic position
        Matches ATLASGAL + SEDIGISM (mm detections) as primary
        Then cross-match with WISE GOs and Hi-GAL for enrichment
        """
        print(f"\n=== Cross-matching catalogs (radius: {match_radius}) ===")

        # Start with primary catalogs that have distances
        primary_catalogs = []
        for name in ['ATLASGAL', 'SEDIGISM']:
            if name in self.catalogs and self.catalogs[name] is not None:
                primary_catalogs.append((name, self.catalogs[name]))

        if not primary_catalogs:
            print("ERROR: No primary catalog (ATLASGAL or SEDIGISM) loaded")
            return None

        # Start with first primary catalog
        name, merged = primary_catalogs[0]
        print(f"Starting with {name}: {len(merged)} sources")

        # Merge additional primary catalogs
        for other_name, other_cat in primary_catalogs[1:]:
            print(f"Merging {other_name}...")
            merged = self._merge_two_catalogs(merged, other_cat, match_radius)

        # Cross-match with WGO
        if 'WGO' in self.catalogs and self.catalogs['WGO'] is not None:
            print(f"Enriching with WISE Green Objects...")
            merged = self._enrich_with_catalog(merged, self.catalogs['WGO'], match_radius, 'WGO_flag')

        # Cross-match with Hi-GAL for luminosity estimates
        if 'HiGAL' in self.catalogs and self.catalogs['HiGAL'] is not None:
            print(f"Enriching with Hi-GAL luminosities...")
            merged = self._enrich_with_catalog(merged, self.catalogs['HiGAL'], match_radius, 'HiGAL_match')

        self.merged = merged
        print(f"\nMerged catalog: {len(merged)} unique sources")
        return merged

    @staticmethod
    def _merge_two_catalogs(cat1, cat2, match_radius):
        """Merge two catalogs, removing duplicates by position"""
        coord1 = SkyCoord(l=cat1['l']*u.deg, b=cat1['b']*u.deg, frame='galactic')
        coord2 = SkyCoord(l=cat2['l']*u.deg, b=cat2['b']*u.deg, frame='galactic')

        idx1, idx2, d2d, d3d = coord1.search_around_sky(coord2, match_radius)

        # Mark matches
        cat1_unmatched = [i for i in range(len(cat1)) if i not in idx1]
        cat2_matched = set(idx2)
        cat2_unmatched = [i for i in range(len(cat2)) if i not in cat2_matched]

        # Combine unmatched sources
        combined = vstack([cat1[cat1_unmatched], cat2[cat2_unmatched]])
        return combined

    @staticmethod
    def _enrich_with_catalog(main_cat, enrich_cat, match_radius, flag_column):
        """Add flag column indicating cross-match with enrichment catalog"""
        coord_main = SkyCoord(l=main_cat['l']*u.deg, b=main_cat['b']*u.deg, frame='galactic')
        coord_enrich = SkyCoord(l=enrich_cat['l']*u.deg, b=enrich_cat['b']*u.deg, frame='galactic')

        idx_main, idx_enrich, d2d, d3d = coord_main.search_around_sky(coord_enrich, match_radius)

        # Add match flag
        main_cat[flag_column] = False
        main_cat[flag_column][idx_main] = True

        # Add luminosity if available and not already present
        if 'L_bol_Lsun' in enrich_cat.colnames and flag_column == 'HiGAL_match':
            if 'L_bol_Lsun_higal' not in main_cat.colnames:
                main_cat['L_bol_Lsun_higal'] = np.nan
                main_cat['L_bol_Lsun_higal'][idx_main] = enrich_cat['L_bol_Lsun'][idx_enrich]

        return main_cat

    def filter_by_distance_luminosity(self, d_min_kpc=0.4, d_max_kpc=2.0, L_min_Lsun=1e4):
        """
        Filter merged catalog by distance and luminosity
        """
        if self.merged is None:
            print("ERROR: Run crossmatch_by_position first")
            return None

        print(f"\n=== Applying filters ===")
        print(f"  Distance: {d_min_kpc} < d < {d_max_kpc} kpc")
        print(f"  Luminosity: L > {L_min_Lsun:.0e} Lsun")

        # Convert distance to kpc if needed
        if 'distance_pc' in self.merged.colnames:
            distance_kpc = self.merged['distance_pc'] / 1000.0
        else:
            print("WARNING: No distance column found")
            return self.merged

        # Apply distance filter
        mask = (distance_kpc >= d_min_kpc) & (distance_kpc <= d_max_kpc)

        # Apply luminosity filter if available
        if 'L_bol_Lsun' in self.merged.colnames:
            valid_L = ~np.isnan(self.merged['L_bol_Lsun'])
            mask &= (valid_L & (self.merged['L_bol_Lsun'] >= L_min_Lsun))
        elif 'L_bol_Lsun_higal' in self.merged.colnames:
            valid_L = ~np.isnan(self.merged['L_bol_Lsun_higal'])
            mask &= (valid_L & (self.merged['L_bol_Lsun_higal'] >= L_min_Lsun))
        else:
            print("WARNING: No luminosity column found; not applying L filter")

        filtered = self.merged[mask]
        print(f"  Retained {len(filtered)} / {len(self.merged)} sources")
        return filtered

    def write_catalog(self, output_filename='myso_candidates_phase1.fits'):
        """Save merged catalog to FITS"""
        output_path = self.data_dir / output_filename
        if self.merged is None:
            print("ERROR: No merged catalog. Run crossmatch_by_position first.")
            return

        self.merged.write(output_path, overwrite=True)
        print(f"Saved catalog to {output_path}")

def main():
    print("\n" + "="*70)
    print("MYSO CATALOG CROSS-MATCH AND FILTERING - PHASE 1")
    print("="*70)

    # Initialize
    merger = MYSOCatalogMerger(DATA_DIR)

    # Load Edenhofer dust map
    dust_data = merger.load_edenhofer_dust_map()

    print("\n" + "="*70)
    print("LOADING CATALOGS")
    print("="*70)

    # Load catalogs
    merger.load_atlasgal()
    merger.load_rms_survey()
    merger.load_sedigism()
    merger.load_wise_green_objects()
    merger.load_higal()

    # Print what we loaded
    print("\n" + "="*70)
    print("CATALOG SUMMARY")
    print("="*70)
    for name, cat in merger.catalogs.items():
        if cat is not None:
            print(f"  ✓ {name:20} {len(cat):6} sources")
        else:
            print(f"  ✗ {name:20} not loaded")

    # Cross-match
    print("\n" + "="*70)
    print("CROSS-MATCHING")
    print("="*70)
    merged = merger.crossmatch_by_position(match_radius=5*u.arcsec)

    if merged is not None:
        print("\n" + "="*70)
        print("APPLYING FILTERS")
        print("="*70)
        # Filter
        filtered = merger.filter_by_distance_luminosity(d_min_kpc=0.4, d_max_kpc=2.0, L_min_Lsun=1e4)

        # Save
        merger.merged = filtered
        merger.write_catalog()

        print("\n" + "="*70)
        print(f"PHASE 1 COMPLETE")
        print("="*70)
        print(f"Final candidate MYSO sample: {len(filtered)} sources")
        if 'catalog_origin' in filtered.colnames:
            print("\nBreakdown by catalog origin:")
            for origin in sorted(set(filtered['catalog_origin'])):
                n = np.sum(filtered['catalog_origin'] == origin)
                print(f"  - {origin:20} {n:6} sources")

        # Distance and luminosity distribution
        if 'distance_pc' in filtered.colnames:
            d_kpc = filtered['distance_pc'] / 1000
            print(f"\nDistance distribution:")
            print(f"  Min: {d_kpc.min():.2f} kpc")
            print(f"  Max: {d_kpc.max():.2f} kpc")
            print(f"  Mean: {np.nanmean(d_kpc):.2f} kpc")

        if 'L_bol_Lsun' in filtered.colnames:
            print(f"\nLuminosity distribution:")
            L_valid = filtered['L_bol_Lsun'][~np.isnan(filtered['L_bol_Lsun'])]
            if len(L_valid) > 0:
                print(f"  N with L estimates: {len(L_valid)}")
                print(f"  Min: {L_valid.min():.2e} Lsun")
                print(f"  Max: {L_valid.max():.2e} Lsun")
                print(f"  Mean: {np.mean(L_valid):.2e} Lsun")

if __name__ == '__main__':
    main()
