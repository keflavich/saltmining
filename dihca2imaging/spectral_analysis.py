#!/usr/bin/env python3
"""
Comprehensive spectral analysis of DIHCA sources following the analysis plan.

For each source in the catalog, this script performs:
1. Create 2-sigma and 10-sigma masks with binary dilation
2. Extract spectra from masked cubes
3. Find significant spectral peaks (5-sigma)
4. Determine source velocities via SIMBAD/line matching
5. Create moment-0 and masked moment images
6. Generate diagnostic image galleries
"""

import os
import glob
import time
import shutil
import numpy as np
import random
from astropy.table import Table, Column
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.stats import mad_std
from astropy.constants import c
from astropy import constants
import astropy.units as u
from astropy import log
from astroquery.simbad import Simbad
from spectral_cube import SpectralCube
from scipy import ndimage
from scipy.ndimage import binary_dilation, label, labeled_comprehension
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)

class SpectralAnalyzer:
    """Main class for analyzing spectral cubes of DIHCA sources"""

    def __init__(self, catalog_file, image_directory, output_directory):
        """
        Initialize the analyzer

        Parameters:
        -----------
        catalog_file : str
            Path to the source catalog FITS file
        image_directory : str
            Directory containing the grouped imaging products
        output_directory : str
            Directory to save analysis results
        """
        self.catalog_file = catalog_file
        self.image_directory = image_directory
        self.output_directory = output_directory

        # Load the source catalog
        self.catalog = Table.read(catalog_file)
        print(f"Loaded catalog with {len(self.catalog)} sources from {len(np.unique(self.catalog['source_field']))} fields")

        # Create output directory structure
        self.setup_output_directories()

        # Initialize SIMBAD for velocity queries
        self.simbad = Simbad()
        self.simbad.add_votable_fields('rvz_radvel')  # radial velocity (correct field name)

        # Common molecular line frequencies (in GHz) for identification
        self.common_lines = {
            'H2CO_3(0,3)-2(0,2)': 218.222,
            'H2CO_3(2,2)-2(2,1)': 218.475,
            'H2CO_3(2,1)-2(2,0)': 218.760,
            'SiO_5-4': 217.105,
            'DCN_3-2': 217.238,
            'C18O_2-1': 219.560,
            '13CO_2-1': 220.399,
            'CH3CN_12-11_0': 220.7472612,
            'CH3CN_12-11_1': 220.743009,
            'CH3CN_12-11_2': 220.73026,
            'CH3CN_12-11_3': 220.709017,
            'SO_6(5)-5(4)': 219.949,
            'H2O_550-643': 232.68670,
            'CH3OH_423-514': 234.68337,
        }

    def setup_output_directories(self):
        """Create the output directory structure"""
        os.makedirs(self.output_directory, exist_ok=True)
        self.spectra_dir = os.path.join(self.output_directory, 'spectra')
        self.moments_dir = os.path.join(self.output_directory, 'moments')
        self.diagnostics_dir = os.path.join(self.output_directory, 'diagnostics')
        self.catalogs_dir = os.path.join(self.output_directory, 'catalogs')

        for directory in [self.spectra_dir, self.moments_dir, self.diagnostics_dir, self.catalogs_dir]:
            os.makedirs(directory, exist_ok=True)

        print(f"Created output directories in {self.output_directory}")

    def load_spectral_cube(self, image_file):
        """
        Load a spectral cube with proper error handling for beam issues

        Parameters:
        -----------
        image_file : str
            Path to the FITS image file

        Returns:
        --------
        cube : SpectralCube
            Loaded spectral cube
        """
        print("  Loading spectral cube...")

        # Load the spectral cube
        # empirically, use_dask is a _lot_ faster for these cubes & operations
        # truncate the edge channels ...
        cube = SpectralCube.read(image_file, use_dask=True)[2:-2]
        print(f"  Cube shape: {cube.shape}")

        # Enable large operations
        cube.allow_huge_operations = True

        return cube

    def find_image_files(self, source_field):
        """Find all spectral cube files for a given source field"""
        pattern = os.path.join(self.image_directory, f"{source_field}_group*_spw*.image.fits")
        files = glob.glob(pattern)
        return sorted(files)

    def create_source_mask(self, cube, source_coord, noise_level):

        # check cube in mask
        ww = cube.wcs.celestial
        ww._naxis = cube.shape[1:]
        if not source_coord.contained_by(ww):
            return

        print(f"Computing cube.max(): ", end="", flush=True)
        t0 = time.time()
        mx = cube.max(axis=0)
        mxmid = np.nanmedian(mx)
        print(f"{time.time()-t0:0.3f} seconds")
        mask_2sig = mx > (2 * noise_level + mxmid)
        mask_1sig = mx > (1 * noise_level + mxmid)

        if mask_1sig.sum() == 0:
            # re-estimate noise
            print(f"Noise level was too high for 1-sigma at {noise_level}; ", end=" ")
            noise_level = np.nanstd(mx)
            print(f"new noise level = {noise_level}")
            mask_2sig = mx > (2 * noise_level + mxmid)
            mask_1sig = mx > (1 * noise_level + mxmid)

        if mask_2sig.sum() == 0:
            print(f"Noise level was too high for 2-sigma at {noise_level}; ", end=" ")
            noise_level = np.nanstd(mx)
            print(f"new noise level = {noise_level}")
            mask_2sig = mx > (2 * noise_level + mxmid)
            mask_1sig = mx > (1 * noise_level + mxmid)
            #raise ValueError(f"CRITICAL ERROR: Found no 2σ pixels in the mask")

        assert mask_2sig.sum() > 0, f"CRITICAL ERROR: Found no 2σ pixels in the mask"
        assert mask_2sig.sum() < mask_1sig.sum()

        dilated_mask = binary_dilation(mask_2sig, mask=mask_1sig, iterations=10)
        if dilated_mask.sum() == 1:
            # weird one-off case, at least so far, where there's a single hot pixel that triggered source creation
            # we can just abandon this "object"
            return

        labeled_mask, num_labels = label(dilated_mask)
        if num_labels > 0:
            source_pixel = cube.wcs.celestial.world_to_pixel(source_coord)
            label_id = labeled_mask[int(source_pixel[1]), int(source_pixel[0])]
            if label_id == 0:
                return
            assert label_id > 0
            source_mask = labeled_mask == label_id
        else:
            raise ValueError(f"CRITICAL ERROR: Found no labeled regions in the mask")

        return source_mask


    def extract_spectrum(self, cube, mask):
        """Extract spectrum from masked cube - MUST succeed or fail loudly
        Refactored to avoid stupid AI wastefulness
        """
        print("    Extracting spectrum from masked cube")
        cube.allow_huge_operations = True # Ensure large operations are enabled

        # Check that mask has valid pixels
        mask_pixels = np.sum(mask)
        if mask_pixels == 0:
            raise ValueError(f"CRITICAL ERROR: Source mask contains 0 pixels! "
                           f"This indicates fundamental problems with:\n"
                           f"  - Data quality (all NaN values?)\n"
                           f"  - Source coordinates (wrong field?)\n"
                           f"  - Noise estimation (too conservative?)\n"
                           f"  - WCS/coordinate transformation\n"
                           f"Mask shape: {mask.shape}")

        # equivalent to 'cube.minimal_subcube(spatial_only=True), but avoids squashing if - as in this case - mask is 2d
        slices = cube.subcube_slices_from_mask(mask, spatial_only=True)

        masked_cube = cube.with_mask(mask)[slices]
        try:
            spectrum = masked_cube.mean(axis=(1, 2), how='slice') # Use how='slice' for memory efficiency
        except TypeError:
            spectrum = masked_cube.mean(axis=(1, 2)) # dask
        print(f"    Extracted spectrum with {len(spectrum)} channels")

        return spectrum

    def find_spectral_peaks(self, spectrum, noise_level, min_sigma=5.0, wing_threshold=1.0):
        """Find spectral peaks - refactored to not do stupid AI loops"""
        print(f"    Finding spectral peaks above {min_sigma}σ")

        flux_data = spectrum

        # Critical check: must have finite data
        finite_mask = np.isfinite(flux_data)
        finite_count = np.sum(finite_mask)

        if finite_count == 0:
            raise ValueError(f"CRITICAL ERROR: Spectrum contains NO finite data for peak finding!\n"
                           f"  - Spectrum length: {len(flux_data)}\n"
                           f"  - All values: {np.unique(flux_data)}\n"
                           f"  - This indicates fundamental data corruption")

        if finite_count < len(flux_data) * 0.1:  # Less than 10% finite data
            raise ValueError(f"CRITICAL ERROR: Spectrum has too few finite values for reliable analysis!\n"
                           f"  - Finite values: {finite_count}/{len(flux_data)} ({100*finite_count/len(flux_data):.1f}%)\n"
                           f"  - Need at least 10% finite data for peak finding")

        print(f"    Spectrum quality check: {finite_count}/{len(flux_data)} finite values")

        # Calculate noise from finite data only
        finite_flux = flux_data[finite_mask]
        if len(finite_flux) < 100:
            print(f"    WARNING: Only {len(finite_flux)} finite values for noise calculation")

        spectrum_noise = mad_std(finite_flux)
        threshold = min_sigma * spectrum_noise

        print(f"    Spectrum noise: {spectrum_noise}")
        print(f"    {min_sigma}σ threshold: {threshold}")
        print(f"    Data range: {np.min(finite_flux)} to {np.max(finite_flux)}")

        # Find peaks above threshold
        peak_mask = spectrum > threshold
        big_mask = spectrum > (wing_threshold * spectrum_noise)
        peak_mask = binary_dilation(peak_mask, mask=big_mask, iterations=10)
        peaks, npeaks = label(peak_mask)

        print(f"    Found {npeaks} peaks above {min_sigma}σ threshold")

        # Create catalog
        peaks_catalog = Table()
        channels = np.arange(len(spectrum))
        lbls = np.arange(1, npeaks + 1)
        if npeaks > 0:
            peaks_catalog['channel'] = labeled_comprehension(channels, peaks, lbls, np.mean, out_dtype=float, default=-999)
            peaks_catalog['peak'] = labeled_comprehension(spectrum, peaks, lbls, np.max, out_dtype=float, default=np.nan)
            peaks_catalog['sum'] = labeled_comprehension(spectrum, peaks, lbls, np.sum, out_dtype=float, default=np.nan)
            peaks_catalog['frequency'] = labeled_comprehension(spectrum.spectral_axis, peaks, lbls, np.mean, out_dtype=float, default=-999)
            peaks_catalog['snr'] = peaks_catalog['peak'] / spectrum_noise

        return peaks_catalog


    def consolidate_peaks(self, peaks_catalog, velocity_threshold=5.0*u.km/u.s):
        """
        Consolidate peaks that are within velocity_threshold of each other.

        Parameters:
        -----------
        peaks_catalog : Table
            Catalog of detected spectral peaks
        velocity_threshold : Quantity
            Velocity threshold (default: 5.0 km/s) for consolidating peaks

        Returns:
        --------
        consolidated_peaks : Table
            Consolidated peaks catalog with unique peaks
        """
        if len(peaks_catalog) == 0:
            return peaks_catalog

        print(f"    Consolidating {len(peaks_catalog)} peaks within {velocity_threshold}")
        velocity_threshold = u.Quantity(velocity_threshold, u.km/u.s)

        # Convert frequencies to velocities if we have frequency information
        if 'frequency' not in peaks_catalog.colnames:
            raise ValueError("No frequency information for velocity consolidation")

        # Sort peaks by frequency for efficient grouping
        # We need to work with relative velocities since rest frequencies are unknown
        sorted_indices = np.argsort(peaks_catalog['frequency'])
        sorted_peaks = peaks_catalog[sorted_indices]

        # Group peaks within velocity_threshold based on relative velocities
        consolidated_peaks = []
        used_indices = set()

        for i in range(len(sorted_peaks)):
            if i in used_indices:
                continue

            # Start a new group with this peak
            group_peaks = [sorted_peaks[i]]
            group_indices = [i]
            ref_freq = sorted_peaks[i]['frequency']

            # Find all peaks within velocity_threshold of this peak
            # Relative velocity: v_rel = (f1 - f2) / f * c
            for j in range(i + 1, len(sorted_peaks)):
                if j in used_indices:
                    continue

                # Calculate relative velocity between the two peaks
                freq_diff = sorted_peaks[j]['frequency'] - ref_freq
                relative_velocity = np.abs((freq_diff / ref_freq) * constants.c).to(u.km/u.s)

                if relative_velocity <= velocity_threshold:
                    group_peaks.append(sorted_peaks[j])
                    group_indices.append(j)
                else:
                    break  # Peaks are sorted by frequency, so no more nearby peaks

            # Mark all peaks in this group as used
            used_indices.update(group_indices)

            # Create consolidated peak (use strongest peak in group)
            snrs = [p['snr'] for p in group_peaks]
            strongest_idx = np.argmax(snrs)
            consolidated_peak = dict(group_peaks[strongest_idx])

            # Add group information
            consolidated_peak['n_peaks_in_group'] = len(group_peaks)
            frequencies = [p['frequency'] for p in group_peaks]
            freq_range = max(frequencies) - min(frequencies)
            # Store frequency range and equivalent velocity range
            consolidated_peak['frequency_range'] = freq_range
            consolidated_peak['velocity_range'] = np.abs((freq_range / ref_freq) * c).to(u.km/u.s)

            consolidated_peaks.append(consolidated_peak)

        # Convert back to Table
        if consolidated_peaks:
            consolidated_catalog = Table(consolidated_peaks)
            print(f"    Consolidated to {len(consolidated_catalog)} unique peaks")

            # Print consolidation summary
            for i, peak in enumerate(consolidated_catalog):
                if peak['n_peaks_in_group'] > 1:
                    v_range = peak['velocity_range']
                    #print(f"      Peak {i+1}: {peak['n_peaks_in_group']} peaks consolidated, "
                    #      f"relative velocity range: {v_range}")
                else:
                    freq = peak['frequency']
                    #print(f"      Peak {i+1}: single peak at f = {freq}")
        else:
            consolidated_catalog = Table()

        return consolidated_catalog

    def query_simbad_velocity(self, source_coord, radius=30*u.arcsec):
        """Query SIMBAD for radial velocity information near the source"""
        print(f"    Querying SIMBAD within {radius} of source")

        # Create a fresh SIMBAD instance to avoid cached field issues
        from astroquery.simbad import Simbad
        simbad = Simbad()
        simbad.add_votable_fields('rvz_radvel')

        result = simbad.query_region(source_coord, radius=radius)

        if result is not None and len(result) > 0:
            print(f"    SIMBAD returned {len(result)} sources")
            print(f"    Available columns: {result.colnames}")

            # Look for radial velocity in various possible column names
            rv_column = None
            for col_name in ['RV_VALUE', 'RVZ_RADVEL', 'rvz_radvel', 'RADVEL', 'RV']:
                if col_name in result.colnames:
                    rv_column = col_name
                    break

            if rv_column:
                print(f"    Using velocity column: {rv_column}")
                rv_mask = ~result[rv_column].mask if hasattr(result[rv_column], 'mask') else ~np.isnan(result[rv_column])
                if np.any(rv_mask):
                    velocities = result[rv_column][rv_mask]
                    mean_velocity = np.mean(velocities)
                    print(f"    Found SIMBAD velocity: {mean_velocity:.1f} km/s from {np.sum(rv_mask)} sources")
                    # SIMBAD velocities are in km/s
                    return u.Quantity(mean_velocity, u.km/u.s)
            else:
                print(f"    No velocity column found in: {result.colnames}")

        print("    No SIMBAD velocity information found")
        return None

    def determine_source_velocity(self, source_coord, peaks_catalog):
        """
        Determine the source velocity using multiple methods

        Parameters:
        -----------
        source_coord : SkyCoord
            Source coordinates
        peaks_catalog : Table
            Catalog of detected spectral peaks

        Returns:
        --------
        velocity : float
            Best estimate of source velocity in km/s
        method : str
            Method used to determine velocity
        """
        print("    Determining source velocity")

        # Method 1: Try SIMBAD
        simbad_velocity = self.query_simbad_velocity(source_coord)
        if simbad_velocity is not None:
            return simbad_velocity, "simbad"

        # Method 2: Try line matching with common molecular lines
        if len(peaks_catalog) > 0 and 'frequency' in peaks_catalog.colnames:
            for line_name, rest_freq in self.common_lines.items():
                for peak in peaks_catalog:
                    # Handle frequency with or without units (convert to GHz)
                    freq_val = peak['frequency']
                    obs_freq_ghz = u.Quantity(freq_val, u.Hz).to(u.GHz)

                    # Calculate velocity using proper units
                    rest_freq_ghz = u.Quantity(rest_freq, u.GHz)
                    velocity = c * (rest_freq_ghz - obs_freq_ghz) / rest_freq_ghz
                    velocity_kms = velocity.to(u.km/u.s)

                    if abs(velocity_kms) < 200 * u.km/u.s:  # Reasonable velocity range
                        print(f"    Matched {line_name}: v = {velocity_kms}")
                        return velocity_kms, f"line_match_{line_name}"

        # Method 3: Use neighbor velocities (placeholder - would need implementation)
        # This would require a more sophisticated approach to find nearby sources

        # Method 4: Default to 0 km/s
        print("    Using default velocity: 0 km/s")
        return 0.0 * u.km/u.s, "default"

    def identify_line_peaks(self, peaks_catalog, source_velocity, velocity_tolerance=20.0*u.km/u.s):
        """
        Identify molecular line peaks after source velocity is determined

        Parameters:
        -----------
        peaks_catalog : Table
            Catalog of detected spectral peaks
        source_velocity : Quantity
            Source velocity (km/s)
        velocity_tolerance : Quantity
            Velocity tolerance for line matching (default: 20 km/s)

        Returns:
        --------
        identified_peaks : Table
            Updated peaks catalog with line identifications added
        """
        print(f"    Identifying line peaks using source velocity: {source_velocity}")

        if len(peaks_catalog) == 0:
            print("    No peaks to identify")
            return peaks_catalog

        if 'frequency' not in peaks_catalog.colnames:
            print("    WARNING: No frequency information for line identification")
            return peaks_catalog

        # Add line_id and line_match_velocity columns if not present
        if 'line_id' not in peaks_catalog.colnames:
            peaks_catalog['line_id'] = [''] * len(peaks_catalog)
        if 'line_match_velocity' not in peaks_catalog.colnames:
            # Store as plain float (km/s assumed)
            peaks_catalog['line_match_velocity'] = np.full(len(peaks_catalog), np.nan)
        if 'line_rest_freq' not in peaks_catalog.colnames:
            # Store as plain float (GHz assumed)
            peaks_catalog['line_rest_freq'] = np.full(len(peaks_catalog), np.nan)

        velocity_tolerance = u.Quantity(velocity_tolerance, u.km/u.s)
        source_velocity = u.Quantity(source_velocity, u.km/u.s)

        n_identified = 0

        # For each peak, try to identify it with known molecular lines
        for i, peak in enumerate(peaks_catalog):
            obs_freq = u.Quantity(peak['frequency'], u.Hz)
            best_match = None
            best_match_velocity_diff = None

            # Try matching with each known line
            for line_name, rest_freq_ghz in self.common_lines.items():
                rest_freq = u.Quantity(rest_freq_ghz, u.GHz)

                # Calculate expected observed frequency if this line is at source velocity
                # obs_freq = rest_freq * (1 - v/c)
                expected_obs_freq = rest_freq * (1 - source_velocity / constants.c)

                # Calculate velocity difference between observed and expected frequency
                velocity_diff = constants.c * (expected_obs_freq - obs_freq) / obs_freq
                velocity_diff = velocity_diff.to(u.km/u.s)

                # Check if within tolerance
                if np.abs(velocity_diff) < velocity_tolerance:
                    # Check if this is the best match so far
                    if best_match is None or np.abs(velocity_diff) < np.abs(best_match_velocity_diff):
                        best_match = line_name
                        best_match_velocity_diff = velocity_diff

            # If we found a match, record it
            if best_match is not None:
                peaks_catalog['line_id'][i] = best_match
                peaks_catalog['line_match_velocity'][i] = best_match_velocity_diff.to_value(u.km/u.s)
                peaks_catalog['line_rest_freq'][i] = u.Quantity(self.common_lines[best_match], u.GHz).to_value(u.GHz)
                n_identified += 1
                print(f"      Peak {i+1} at {obs_freq.to(u.GHz):.4f}: Identified as {best_match} "
                      f"(Δv = {best_match_velocity_diff:.2f})")
            else:
                peaks_catalog['line_id'][i] = 'unidentified'
                print(f"      Peak {i+1} at {obs_freq.to(u.GHz):.4f}: No line identification")

        print(f"    Identified {n_identified}/{len(peaks_catalog)} line peaks")

        return peaks_catalog

    def create_moment_images_for_peaks(self, cube, mask, consolidated_peaks, line_width=5.0*u.km/u.s):
        """
        Create moment-0 and moment-1 images for each consolidated peak

        Parameters:
        -----------
        cube : SpectralCube
            The spectral cube
        mask : numpy.ndarray
            3D mask for the source
        consolidated_peaks : Table
            Consolidated peaks catalog
        line_width : Quantity
            Line width for moment calculation (default: 5.0 km/s)

        Returns:
        --------
        all_moments : dict
            Dictionary containing moment images for each peak
        """
        print(f"    Creating moment images for {len(consolidated_peaks)} consolidated peaks")

        all_moments = {}

        if len(consolidated_peaks) == 0:
            print("    No peaks found - creating default moment images from the middle frequency of the cube")
            # Create default moment images using v=0, but return in consistent nested structure
            default_moments = self.create_moment_images(cube, mask, velocity=0.0 * u.km/u.s, rest_freq=cube.spectral_axis.mean(), line_width=line_width)
            # Return empty dict instead to be consistent with the nested structure expectation
            return {}

        for i, peak in enumerate(consolidated_peaks):
            peak_id = f"peak_{i+1}"
            print(f"    Creating moments for {peak_id} at v = {peak['velocity'] if 'velocity' in peak else 0.0 * u.km/u.s}")

            # Get frequency (preserve units, convert to Hz)
            freq_val = peak['frequency']
            rest_freq = u.Quantity(freq_val, u.Hz)

            # Create moment images for this peak
            moments = self.create_moment_images(cube, mask,
                                              velocity=peak['velocity'] if 'velocity' in peak else 0.0 * u.km/u.s,
                                              rest_freq=rest_freq,
                                              line_width=line_width)

            # Store with peak identifier (keep units)
            all_moments[peak_id] = {
                'moments': moments,
                'velocity': peak['velocity'] if 'velocity' in peak else 0.0 * u.km/u.s,
                'frequency': rest_freq,  # Now has units
                'snr': peak['snr'],
                'n_peaks_in_group': peak.get('n_peaks_in_group', 1)
            }

        return all_moments

    def create_moment_images(self, cube, mask, velocity, rest_freq, line_width=5.0*u.km/u.s):
        """
        Create moment-0 and moment-1 images

        Parameters:
        -----------
        cube : SpectralCube
            The spectral cube
        mask : numpy.ndarray
            3D mask for the source
        velocity : Quantity
            Source velocity (with units of km/s or equivalent)
        rest_freq : Quantity
            Rest frequency (with units of Hz or equivalent)
        line_width : Quantity
            Line width for moment calculation (default: 5.0 km/s)

        Returns:
        --------
        moments : dict
            Dictionary containing moment images
        """
        print(f"    Creating moment images around freq = {rest_freq} and v = {velocity}")
        t0 = time.time()
        assert rest_freq.unit.is_equivalent(u.Hz)

        moments = {}

        # Convert velocity to km/s if needed
        v_center = u.Quantity(velocity, u.km/u.s)
        v_range = u.Quantity(line_width, u.km/u.s)

        # Convert to frequency
        freq_center = rest_freq * (1 - v_center / constants.c)
        freq_range = rest_freq * v_range / constants.c

        # Select velocity range
        freq_min = freq_center - freq_range
        freq_max = freq_center + freq_range

        # Create subcube
        slices = cube.subcube_slices_from_mask(mask, spatial_only=True)
        masked_subcube = cube.with_mask(mask).spectral_slab(freq_min, freq_max)[slices]

        # Create moments
        moments['moment0'] = masked_subcube.moment0()
        moments['moment1'] = masked_subcube.moment1()
        moments['peak_intensity'] = masked_subcube.max(axis=0)
        try:
            moments['velocity_of_peak'] = masked_subcube.argmax_world(axis=0)
        except IndexError:
            # Hack around https://github.com/radio-astro-tools/spectral-cube/issues/982.  Not good, but not expected to be a common case (only affects size-1 slices)
            moments['velocity_of_peak'] = moments['moment1']

        print(f"    Created moments from {masked_subcube.shape[0]} channels in image size {masked_subcube.shape[1:]}.  Time={time.time()-t0:0.3f} seconds")

        return moments

    def create_masked_moment_images_for_peaks(self, cube, spatial_mask, consolidated_peaks, line_width=5.0*u.km/u.s):
        """
        Create masked moment images for each consolidated peak using adaptive thresholding

        Parameters:
        -----------
        cube : SpectralCube
            The spectral cube
        consolidated_peaks : Table
            Consolidated peaks catalog
        line_width : Quantity
            Line width for moment calculation (default: 5.0 km/s)

        Returns:
        --------
        all_masked_moments : dict
            Dictionary containing masked moment images for each peak
        """
        print(f"    Creating masked moment images for {len(consolidated_peaks)} consolidated peaks")

        all_masked_moments = {}

        if len(consolidated_peaks) == 0:
            print("    No peaks found - creating default masked moment images from the middle frequency of the cube")
            # Return empty dict instead to be consistent with the nested structure expectation
            return {}

        for i, peak in enumerate(consolidated_peaks):
            peak_id = f"peak_{i+1}"
            print(f"    Creating masked moments for {peak_id} at v = {peak['velocity'] if 'velocity' in peak else 0.0 * u.km/u.s}")

            # Get frequency (preserve units, convert to Hz)
            freq_val = peak['frequency']
            rest_freq = u.Quantity(freq_val, u.Hz)

            # Create masked moment images for this peak
            masked_moments = self.create_masked_moment_images(cube,
                                                              spatial_mask=spatial_mask,
                                                              velocity=peak['velocity'] if 'velocity' in peak else 0.0 * u.km/u.s,
                                                              rest_freq=rest_freq,
                                                              line_width=line_width)

            # Store with peak identifier (keep units)
            all_masked_moments[peak_id] = {
                'masked_moments': masked_moments,
                'velocity': peak['velocity'] if 'velocity' in peak else 0.0 * u.km/u.s,
                'frequency': rest_freq,
                'snr': peak['snr'],
                'n_peaks_in_group': peak.get('n_peaks_in_group', 1)
            }

        return all_masked_moments

    def create_masked_moment_images(self, cube, spatial_mask, velocity, rest_freq, line_width=5.0*u.km/u.s):
        """
        Create masked moment images using adaptive thresholding

        Parameters:
        -----------
        cube : SpectralCube
            The spectral cube
        velocity : Quantity
            Source velocity (with units of km/s or equivalent)
        rest_freq : Quantity
            Rest frequency (with units of Hz or equivalent)
        line_width : Quantity
            Line width for moment calculation (default: 5.0 km/s)

        Returns:
        --------
        masked_moments : dict
            Dictionary containing masked moment images
        """
        print(f"    Creating masked moment images around freq = {rest_freq} and v = {velocity}. ", end="", flush=True)
        t0 = time.time()

        masked_moments = {}

        # Similar to create_moment_images but with adaptive masking
        v_center = u.Quantity(velocity, u.km/u.s)
        v_range = 2 * u.Quantity(line_width, u.km/u.s)  # ±2 FWHM

        freq_center = v_center.to(u.GHz, u.doppler_radio(rest_freq))
        freq_range = rest_freq * v_range / constants.c

        freq_min = freq_center - freq_range
        freq_max = freq_center + freq_range

        # be a little generous with the assertions around the edges
        assert freq_min > cube.spectral_axis.min() - freq_range*2
        assert freq_max < cube.spectral_axis.max() + freq_range*2

        slices = cube.subcube_slices_from_mask(spatial_mask, spatial_only=True)
        subcube = cube.spectral_slab(freq_min, freq_max).with_mask(spatial_mask)[slices]

        noise_level = subcube.mad_std()

        # Create 1-sigma and 5-sigma masks, with iterative reduction if needed
        for sigma_factor in [5.0, 4.0, 3.0, 2.0, 1.5]:
            threshold = sigma_factor * noise_level

            # Use spectral-cube comparison directly
            high_sigma_mask = subcube > threshold

            if np.any(high_sigma_mask.include()):
                npix = int(np.sum(high_sigma_mask.include()))
                print(f"    Using {sigma_factor}σ threshold: {npix} pixels")
                break
        else:
            raise ValueError("No significant pixels found above 1.5 sigma")

        # Create 1-sigma mask using spectral-cube comparison
        low_sigma_mask = subcube > noise_level

        # Grow high-sigma mask into low-sigma mask
        if np.any(high_sigma_mask.include()):
            # Convert spectral-cube masks to numpy arrays for dilation operations
            high_sigma_array = high_sigma_mask.include()
            low_sigma_array = low_sigma_mask.include()

            final_mask = binary_dilation(high_sigma_array, mask=low_sigma_array, iterations=10)
        else:
            peak_ind = subcube.argmax()
            final_mask = low_sigma_mask.include()*0
            low_sigma_array = low_sigma_mask.include()
            final_mask[peak_ind] = 1
            final_mask = binary_dilation(final_mask, mask=low_sigma_array, iterations=10)

        # Apply mask and create moments
        # Convert numpy mask to proper spectral-cube mask
        from spectral_cube import BooleanArrayMask

        # Create a proper spectral-cube mask
        mask_obj = BooleanArrayMask(final_mask, wcs=subcube.wcs)
        masked_subcube = subcube.with_mask(mask_obj)

        masked_moments['masked_moment0'] = masked_subcube.moment0()
        masked_moments['masked_moment1'] = masked_subcube.moment1()
        try:
            masked_moments['velocity_of_peak'] = masked_subcube.argmax_world(axis=0)
        except IndexError:
            # Hack around https://github.com/radio-astro-tools/spectral-cube/issues/982.  Not good, but not expected to be a common case (only affects size-1 slices)
            masked_moments['velocity_of_peak'] = masked_moments['moment1']

        print(f"    Created masked moments with {np.sum(final_mask)} masked pixels in t={time.time()-t0:0.3f} seconds")

        return masked_moments

    def create_diagnostic_gallery(self, source_id, spectrum, moments, masked_moments, peaks_catalog, output_path, peak_id=None, peak_info=None):
        """
        Create a diagnostic image gallery for a specific spectral line/peak

        Parameters:
        -----------
        source_id : str
            Source identifier
        spectrum : object
            Extracted spectrum
        moments : dict
            Dictionary of moment images for this specific peak
        masked_moments : dict
            Dictionary of masked moment images for this specific peak
        peaks_catalog : Table
            Catalog entry for this specific peak (single-row table)
        output_path : str
            Path to save the gallery image
        peak_id : str, optional
            Peak identifier (e.g., "peak_1", "peak_2")
        peak_info : dict, optional
            Additional information about this peak (velocity, frequency, SNR, etc.)
        """
        print(f"    Creating diagnostic gallery: {output_path}")
        if peak_id:
            print(f"    Peak ID: {peak_id}")
            if peak_info:
                print(f"    Peak info: velocity={peak_info.get('velocity', 'N/A')}, "
                      f"frequency={peak_info.get('frequency', 'N/A')}, SNR={peak_info.get('snr', 'N/A')}")
        print(f"    Spectrum from {spectrum.spectral_axis.min()} to {spectrum.spectral_axis.max()}")

        try:
            root_source_id = "_".join(source_id.split("_")[:-1])
            continuum_image_name = f'/orange/adamginsburg/salt/dihca/{root_source_id}.quick_cont.pbclean.image.fits'
            conthdu = fits.open(continuum_image_name)[0]
            print(f"Found {continuum_image_name}")
        except Exception as ex:
            print(ex)

        fig = plt.figure(figsize=(15, 12))
        gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)

        # Row 1: moment-0, peak-intensity, masked moment-0
        ax1 = fig.add_subplot(gs[0, 0])
        if moments and moments.get('moment0') is not None:
            moment0_data = moments['moment0'].value if hasattr(moments['moment0'], 'value') else moments['moment0']
            if np.any(np.isfinite(moment0_data)):
                im1 = ax1.imshow(moment0_data, origin='lower', cmap='viridis')
                plt.colorbar(im1, ax=ax1, shrink=0.8)
                ax1.set_title('Moment-0')
            else:
                ax1.text(0.5, 0.5, 'No finite data', ha='center', va='center', transform=ax1.transAxes)
                ax1.set_title('Moment-0 (no data)')
        else:
            ax1.text(0.5, 0.5, 'No moment-0 data', ha='center', va='center', transform=ax1.transAxes)
            ax1.set_title('Moment-0 (unavailable)')

        ax2 = fig.add_subplot(gs[0, 1])
        if moments and moments.get('peak_intensity') is not None:
            peak_data = moments['peak_intensity'].value if hasattr(moments['peak_intensity'], 'value') else moments['peak_intensity']
            if np.any(np.isfinite(peak_data)):
                im2 = ax2.imshow(peak_data, origin='lower', cmap='plasma')
                plt.colorbar(im2, ax=ax2, shrink=0.8)
                ax2.set_title('Peak Intensity')
            else:
                ax2.text(0.5, 0.5, 'No finite data', ha='center', va='center', transform=ax2.transAxes)
                ax2.set_title('Peak Intensity (no data)')
        else:
            ax2.text(0.5, 0.5, 'No peak intensity data', ha='center', va='center', transform=ax2.transAxes)
            ax2.set_title('Peak Intensity (unavailable)')

        ax3 = fig.add_subplot(gs[0, 2])
        if masked_moments and masked_moments.get('masked_moment0') is not None:
            masked_data = masked_moments['masked_moment0'].value if hasattr(masked_moments['masked_moment0'], 'value') else masked_moments['masked_moment0']
            if np.any(np.isfinite(masked_data)):
                im3 = ax3.imshow(masked_data, origin='lower', cmap='viridis')
                plt.colorbar(im3, ax=ax3, shrink=0.8)
                ax3.set_title('Masked Moment-0')
            else:
                ax3.text(0.5, 0.5, 'No finite data', ha='center', va='center', transform=ax3.transAxes)
                ax3.set_title('Masked Moment-0 (no data)')
        else:
            ax3.text(0.5, 0.5, 'No masked moment-0 data', ha='center', va='center', transform=ax3.transAxes)
            ax3.set_title('Masked Moment-0 (unavailable)')

        # Row 2: moment-1, velocity-of-peak, masked moment-1
        ax4 = fig.add_subplot(gs[1, 0])
        if moments and moments.get('moment1') is not None:
            moment1_obj = moments['moment1']

            # Convert moment-1 from frequency to velocity if peak_info is available
            if peak_info and 'frequency' in peak_info and hasattr(moment1_obj, 'unit'):
                rest_freq = peak_info['frequency']
                moment1_velocity = moment1_obj.to(u.km/u.s, u.doppler_radio(rest_freq))

                moment1_display = moment1_velocity
                moment1_label = 'Moment-1'
                moment1_unit_label = 'km/s'
            else:
                moment1_display = moment1_obj
                moment1_label = 'Moment-1'
                moment1_unit_label = str(moment1_obj.unit) if hasattr(moment1_obj, 'unit') else ''

            # Extract data for display (matplotlib needs arrays, not Quantity objects)
            display_data = moment1_display.value if hasattr(moment1_display, 'value') else moment1_display

            if np.any(np.isfinite(display_data)):
                im4 = ax4.imshow(display_data, origin='lower', cmap='RdBu_r')
                cbar4 = plt.colorbar(im4, ax=ax4, shrink=0.8)
                if moment1_unit_label:
                    cbar4.set_label(moment1_unit_label)
                ax4.set_title(moment1_label)
            else:
                ax4.text(0.5, 0.5, 'No finite data', ha='center', va='center', transform=ax4.transAxes)
                ax4.set_title('Moment-1 (no data)')
        else:
            ax4.text(0.5, 0.5, 'No moment-1 data', ha='center', va='center', transform=ax4.transAxes)
            ax4.set_title('Moment-1 (unavailable)')

        ax5 = fig.add_subplot(gs[1, 1])
        ax5.text(0.5, 0.5, 'Velocity of Peak', ha='center', va='center', transform=ax5.transAxes)
        ax5.set_title('Velocity of Peak')

        ax6 = fig.add_subplot(gs[1, 2])
        if masked_moments and masked_moments.get('masked_moment1') is not None:
            masked_moment1_obj = masked_moments['masked_moment1']

            # Convert masked moment-1 from frequency to velocity if peak_info is available
            if peak_info and 'frequency' in peak_info and hasattr(masked_moment1_obj, 'unit'):
                rest_freq = peak_info['frequency']
                masked_moment1_velocity = masked_moment1_obj.to(u.km/u.s, u.doppler_radio(rest_freq))

                masked_moment1_display = masked_moment1_velocity
                masked_moment1_label = 'Masked Moment-1'
                masked_moment1_unit_label = 'km/s'
            else:
                masked_moment1_display = masked_moment1_obj
                masked_moment1_label = 'Masked Moment-1'
                masked_moment1_unit_label = str(masked_moment1_obj.unit) if hasattr(masked_moment1_obj, 'unit') else ''

            # Extract data for display (matplotlib needs arrays, not Quantity objects)
            display_data = masked_moment1_display.value if hasattr(masked_moment1_display, 'value') else masked_moment1_display

            if np.any(np.isfinite(display_data)):
                im6 = ax6.imshow(display_data, origin='lower', cmap='RdBu_r')
                cbar6 = plt.colorbar(im6, ax=ax6, shrink=0.8)
                if masked_moment1_unit_label:
                    cbar6.set_label(masked_moment1_unit_label)
                ax6.set_title(masked_moment1_label)
            else:
                ax6.text(0.5, 0.5, 'No finite data', ha='center', va='center', transform=ax6.transAxes)
                ax6.set_title('Masked Moment-1 (no data)')
        else:
            ax6.text(0.5, 0.5, 'No masked moment-1 data', ha='center', va='center', transform=ax6.transAxes)
            ax6.set_title('Masked Moment-1 (unavailable)')

        # Row 3: spectrum (spans all columns)
        ax7 = fig.add_subplot(gs[2, :])
        if spectrum is not None:
            # SpectralCube spectrum object
            x_axis = spectrum.spectral_axis
            y_axis = spectrum.quantity # drop the spectral-cube metadata, it cracks something at the mpl level
            try:
                ax7.plot(x_axis, y_axis, 'k-', alpha=0.7)
                ax7.set_xlabel('Frequency/Velocity')
                ax7.set_ylabel('Flux')
            except ValueError as ex:
                print("MAJOR FAILURE that we're skipping: can't plot a simple spectrum")

            # Mark the current peak being displayed with emphasis
            if len(peaks_catalog) > 0:
                for peak in peaks_catalog:
                    if 'frequency' in peaks_catalog.colnames:
                        peak_x = peak['frequency']
                    else:
                        raise ValueError(f"No frequency column in peaks catalog")
                    # Highlight this peak more prominently since it's the one being displayed
                    ax7.axvline(peak_x, color='red', alpha=0.9, linestyle='-', linewidth=2, label=peak_x)

                    # Build label with SNR and line ID if available
                    label = f"SNR={peak['snr']:.1f}"
                    if 'line_id' in peaks_catalog.colnames and peak['line_id']:
                        line_id = peak['line_id']
                        if line_id != 'unidentified':
                            label = f"{line_id} {peak_x}\n{label}"

                    ax7.text(peak_x, peak['peak_intensity'] if 'peak_intensity' in peak else peak['peak'], label,
                            rotation=90, ha='right', va='bottom', fontsize=10,
                            bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))
        else:
            ax7.text(0.5, 0.5, 'No spectrum data available', ha='center', va='center', transform=ax7.transAxes)
            ax7.set_title('Extracted Spectrum (unavailable)')

        # Create title with peak-specific information
        title = f'Diagnostic Gallery: {source_id}'
        if peak_id:
            title += f' - {peak_id}'
            if peak_info:
                # Add line identification if available
                if 'line_id' in peaks_catalog.colnames and len(peaks_catalog) > 0:
                    line_id = peaks_catalog[0]['line_id']
                    if line_id and line_id != 'unidentified':
                        title += f' ({line_id})'
                # Add velocity and SNR info
                velocity = peak_info.get('velocity')
                snr = peak_info.get('snr')
                if velocity is not None and snr is not None:
                    vel_val = velocity.to(u.km/u.s).value if hasattr(velocity, 'to') else velocity
                    title += f' | v={vel_val:.1f} km/s, SNR={snr:.1f}'

        plt.suptitle(title, fontsize=16)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"    Saved diagnostic gallery successfully")


    def analyze_source(self, source_row, image_files=None):
        """
        Analyze a single source following the complete analysis plan

        Parameters:
        -----------
        source_row : astropy.table.Row
            Row from the source catalog containing source information
        """
        source_id = source_row['source_id']
        source_field = source_row['source_field']

        print(f"\n=== Analyzing {source_id} ===")

        # Create source coordinate
        source_coord = SkyCoord(ra=source_row['ra_deg']*u.deg,
                               dec=source_row['dec_deg']*u.deg,
                               frame='icrs')

        # Find image files for this source field (or use provided subset)
        if image_files is None:
            image_files = self.find_image_files(source_field)
        if not image_files:
            print(f"No image files found for {source_field}")
            return

        print(f"Found {len(image_files)} image files for {source_field}")

        # Create source-specific output directory
        source_output_dir = os.path.join(self.output_directory, source_id)
        os.makedirs(source_output_dir, exist_ok=True)

        # Initialize results storage
        all_results = []

        # Process each spectral window
        for image_file in image_files:
            print(f"\n--- Processing {os.path.basename(image_file)} ---")

            # Load the spectral cube
            cube = self.load_spectral_cube(image_file)

            # CRITICAL: Check if we have any finite data at all
            spw_name = os.path.basename(image_file).replace('.image.fits', '')
            print(f"  Validating cube data quality for {spw_name}...")
            # Use spectral-cube methods properly
            # sample subcube has to be cut out of the middle, as there are edges that are junk but there can't be junk centers
            sample_subcube = cube[cube.shape[0]//2-50:cube.shape[0]//2+50, cube.shape[1]//2-50:cube.shape[1]//2+50, cube.shape[2]//2-50:cube.shape[2]//2+50]

            mx = sample_subcube.max()
            if mx == 0*cube.unit or np.isnan(mx):
                cubemx = cube.max()
                if np.isnan(cubemx):
                    print(f"  CRITICAL ERROR: Spectral cube contains only NaN/inf values - deleting file {image_file}")
                    if image_file.endswith(".fits"):
                        os.remove(image_file)
                        image_file = image_file.replace(".fits", "")
                        shutil.rmtree(image_file)
                    else:
                        shutil.rmtree(image_file)
                        os.remove(image_file + ".fits")
                    continue
                else:
                    raise ValueError(f"CRITICAL ERROR: Spectral cube contains only NaN/inf values in sample but is not all nan")

            # Calculate noise from actual finite data
            noise_level = sample_subcube.mad_std()
            print(f"  Noise from finite data: {noise_level}")

            # Step 1: Create source mask
            mask = self.create_source_mask(cube, source_coord, noise_level)
            if mask is None:
                print(f"Source {source_id} is not detected in {os.path.basename(image_file)}")
                gallery_path = os.path.join(source_output_dir, f"{source_id}_{spw_name}_no_peaks_diagnostics.png")
                print(f"Computing max spectrum ", end="", flush=True)
                t0 = time.time()
                spectrum = cube.max(axis=(1,2))
                print(f" in t={time.time() - t0:0.1f} seconds")
                self.create_diagnostic_gallery(source_id, spectrum, {}, {},
                                               {}, gallery_path)
                continue
            else:
                if 'detections' in source_row and isinstance(source_row['detections'], list):
                    source_row['detections'].append(os.path.basename(image_file))
                else:
                    source_row['detections'] = [os.path.basename(image_file)]

            # Step 2: Extract spectrum
            spectrum = self.extract_spectrum(cube, mask)
            if np.all(np.isnan(spectrum)):
                raise ValueError(f"CRITICAL ERROR: Extracted spectrum contains only NaN values!\n"
                               f"  - Spectrum length: {len(spectrum)}\n"
                               f"  - This indicates the spectral cube data is fundamentally broken")

            # Step 3: Find spectral peaks
            peaks_catalog = self.find_spectral_peaks(spectrum, noise_level)

            # Step 3.5: Consolidate nearby peaks within 5 km/s
            consolidated_peaks = self.consolidate_peaks(peaks_catalog, velocity_threshold=5.0*u.km/u.s)

            # Step 4: Determine source velocity (using consolidated peaks)
            velocity, velocity_method = self.determine_source_velocity(source_coord, consolidated_peaks)

            # Step 4.5: Attempt to identify the remaining line peaks
            consolidated_peaks = self.identify_line_peaks(consolidated_peaks, velocity)

            # Step 5: Create moment images for each consolidated peak
            all_moments = self.create_moment_images_for_peaks(cube, mask, consolidated_peaks)

            # Step 6: Create masked moment images for each consolidated peak
            all_masked_moments = self.create_masked_moment_images_for_peaks(cube, mask, consolidated_peaks)

            # Step 7: Create diagnostic galleries - one for each detected peak

            # Get all unique peak_ids from both all_moments and all_masked_moments
            peak_ids = set()
            if all_moments:
                peak_ids.update(all_moments.keys())
            if all_masked_moments:
                peak_ids.update(all_masked_moments.keys())

            # Create a diagnostic gallery for each detected peak
            if peak_ids:
                for peak_id in sorted(peak_ids):
                    # Get moments for this specific peak
                    if peak_id in all_moments:
                        display_moments = all_moments[peak_id]['moments']
                        peak_info = all_moments[peak_id]
                    else:
                        display_moments = {}
                        peak_info = {}

                    # Get masked moments for this specific peak
                    if peak_id in all_masked_moments and all_masked_moments[peak_id].get('masked_moments'):
                        display_masked_moments = all_masked_moments[peak_id]['masked_moments']
                    else:
                        display_masked_moments = {}

                    # Get the corresponding entry from consolidated_peaks for this peak
                    peak_idx = int(peak_id.split('_')[1]) - 1  # Extract peak number (peak_1 -> index 0)
                    if 0 <= peak_idx < len(consolidated_peaks):
                        peak_catalog_entry = consolidated_peaks[peak_idx:peak_idx+1]  # Single-row table
                    else:
                        peak_catalog_entry = Table()  # Empty table

                    # Create gallery filename with peak identifier
                    gallery_path = os.path.join(source_output_dir,
                                              f"{source_id}_{spw_name}_{peak_id}_diagnostics.png")

                    # Create gallery for this peak with its specific catalog entry
                    self.create_diagnostic_gallery(source_id, spectrum, display_moments,
                                                 display_masked_moments, peak_catalog_entry,
                                                 gallery_path, peak_id=peak_id, peak_info=peak_info)
                    print(f"    Created diagnostic gallery for {peak_id}: {os.path.basename(gallery_path)}")
            else:
                # No peaks found - create a single default gallery
                gallery_path = os.path.join(source_output_dir, f"{source_id}_{spw_name}_no_peaks_diagnostics.png")
                self.create_diagnostic_gallery(source_id, spectrum, {}, {},
                                             consolidated_peaks, gallery_path)

            # Save spectrum first
            spectrum_file = os.path.join(source_output_dir, f"{source_id}_{spw_name}_spectrum.fits")
            if hasattr(spectrum, 'write'):
                spectrum.write(spectrum_file, overwrite=True)
                print(f"    Saved spectrum: {os.path.basename(spectrum_file)}")

            # Create one result row for each peak
            # Get all unique peak_ids from both all_moments and all_masked_moments
            peak_ids = set()
            if all_moments:
                peak_ids.update(all_moments.keys())
            if all_masked_moments:
                peak_ids.update(all_masked_moments.keys())

            # If no peaks found, create a single result row with no peak info
            if len(peak_ids) == 0:
                result = {
                    'source_id': source_id,
                    'spw': spw_name,
                    'image_file': image_file,
                    'spectrum_file': spectrum_file,
                    'peak_id': 'none',
                    'peak_velocity': 0.0 * u.km/u.s,
                    'peak_frequency': 0.0 * u.Hz,
                    'peak_snr': 0.0,
                    'n_peaks_in_group': 0,
                    'n_original_peaks': len(peaks_catalog),
                    'n_consolidated_peaks': 0,
                    'noise_level': noise_level,
                    'velocity_method': velocity_method,
                    'line_id': '',
                    'line_match_velocity': np.nan * u.km/u.s,
                    'line_rest_freq': np.nan * u.GHz
                }
                all_results.append(result)

            # Create one row for each peak
            for peak_id in sorted(peak_ids):
                result = {
                    'source_id': source_id,
                    'spw': spw_name,
                    'image_file': image_file,
                    'spectrum_file': spectrum_file,
                    'peak_id': peak_id,
                    'n_original_peaks': len(peaks_catalog),
                    'n_consolidated_peaks': len(consolidated_peaks),
                    'noise_level': noise_level,
                    'velocity_method': velocity_method
                }

                # Add moment information for this peak
                if all_moments and peak_id in all_moments:
                    peak_data = all_moments[peak_id]
                    result['peak_velocity'] = peak_data.get('velocity', 0.0 * u.km/u.s)
                    result['peak_frequency'] = peak_data.get('frequency', 0.0 * u.Hz)
                    result['peak_snr'] = peak_data.get('snr', 0.0)
                    result['n_peaks_in_group'] = peak_data.get('n_peaks_in_group', 1)

                # Add line identification information from consolidated_peaks
                peak_idx = int(peak_id.split('_')[1]) - 1  # Extract peak number
                if 0 <= peak_idx < len(consolidated_peaks):
                    peak_row = consolidated_peaks[peak_idx]
                    result['line_id'] = peak_row.get('line_id', '')
                    # line_match_velocity and line_rest_freq are stored as plain floats in the table
                    result['line_match_velocity'] = peak_row.get('line_match_velocity', np.nan) * u.km/u.s
                    result['line_rest_freq'] = peak_row.get('line_rest_freq', np.nan) * u.GHz
                else:
                    result['line_id'] = ''
                    result['line_match_velocity'] = np.nan * u.km/u.s
                    result['line_rest_freq'] = np.nan * u.GHz

                    # Save moment images for this peak
                    if peak_data.get('moments'):
                        moments = peak_data['moments']
                        peak_velocity = peak_data['velocity']
                        # Extract value for filename (acceptable use of .value for formatting)
                        vel_val = peak_velocity.to(u.km/u.s).value if hasattr(peak_velocity, 'to') else peak_velocity

                        for moment_type in ['moment0', 'moment1']:
                            if moment_type in moments and moments[moment_type] is not None:
                                moment_file = os.path.join(source_output_dir,
                                                         f"{source_id}_{spw_name}_{peak_id}_{moment_type}_v{vel_val:.1f}.fits")
                                moments[moment_type].write(moment_file, overwrite=True)
                                result[f'{moment_type}_file'] = moment_file
                                print(f"    Saved {moment_type} for {peak_id}: {os.path.basename(moment_file)}")

                # Add masked moment information for this peak
                if all_masked_moments and peak_id in all_masked_moments:
                    peak_data = all_masked_moments[peak_id]

                    # If we didn't get velocity from all_moments, get it from all_masked_moments
                    if 'peak_velocity' not in result:
                        result['peak_velocity'] = peak_data.get('velocity', 0.0 * u.km/u.s)
                        result['peak_frequency'] = peak_data.get('frequency', 0.0 * u.Hz)
                        result['peak_snr'] = peak_data.get('snr', 0.0)
                        result['n_peaks_in_group'] = peak_data.get('n_peaks_in_group', 1)

                    # Save masked moment images for this peak
                    if peak_data.get('masked_moments'):
                        masked_moments = peak_data['masked_moments']
                        peak_velocity = peak_data['velocity']
                        # Extract value for filename (acceptable use of .value for formatting)
                        vel_val = peak_velocity.to(u.km/u.s).value if hasattr(peak_velocity, 'to') else peak_velocity

                        for moment_type in ['masked_moment0', 'masked_moment1', 'velocity_of_peak']:
                            if moment_type in masked_moments and masked_moments[moment_type] is not None:
                                moment_file = os.path.join(source_output_dir,
                                                         f"{source_id}_{spw_name}_{peak_id}_{moment_type}_v{vel_val:.1f}.fits")
                                masked_moments[moment_type].write(moment_file, overwrite=True)
                                result[f'{moment_type}_file'] = moment_file
                                print(f"    Saved {moment_type} for {peak_id}: {os.path.basename(moment_file)}")

                all_results.append(result)

            # Save consolidated peaks catalog
            if len(consolidated_peaks) > 0:
                peaks_file = os.path.join(source_output_dir, f"{source_id}_{spw_name}_consolidated_peaks.fits")
                consolidated_peaks.write(peaks_file, overwrite=True)
                result['peaks_file'] = peaks_file
                print(f"    Saved consolidated peaks catalog: {os.path.basename(peaks_file)}")

            # Also save original peaks catalog for comparison
            if len(peaks_catalog) > 0:
                orig_peaks_file = os.path.join(source_output_dir, f"{source_id}_{spw_name}_original_peaks.fits")
                peaks_catalog.write(orig_peaks_file, overwrite=True)
                result['original_peaks_file'] = orig_peaks_file
                print(f"    Saved original peaks catalog: {os.path.basename(orig_peaks_file)}")
            print(f"  Completed analysis of {spw_name}")

        # Save summary results for this source
        if all_results:
            summary_table = Table(all_results)
            summary_file = os.path.join(source_output_dir, f"{source_id}_analysis_summary.fits")
            summary_table.write(summary_file, overwrite=True)
            print(f"Saved analysis summary to {summary_file}")

        print(f"=== Completed analysis of {source_id} ===")

    def analyze_all_sources(self, max_sources=None):
        """
        Analyze all sources in the catalog

        Parameters:
        -----------
        max_sources : int, optional
            Maximum number of sources to analyze (for testing)
        """
        print(f"Starting analysis of {len(self.catalog)} sources...")

        self.catalog['detections'] = None
        sources_to_analyze = self.catalog[:max_sources] if max_sources else self.catalog
        # shuffle is just because I've done the first half a billion times now...
        indices = np.arange(len(sources_to_analyze))
        np.random.shuffle(indices)

        print(f"Actually, starting analysis of {len(sources_to_analyze)} sources...")

        for i, source_row in enumerate(sources_to_analyze[indices]):
            print(f"\n{'='*60}")
            print(f"Progress: {i+1}/{len(sources_to_analyze)}")
            print(f"{'='*60}")

            t0 = time.time()
            # Analyze this source - all errors will raise and halt execution
            self.analyze_source(source_row)
            print(f"Completed source {source_row['source_field']} {source_row['source_number']} in t={time.time()-t0:0.1f} seconds")

        print(f"\n{'='*60}")
        print("Analysis complete!")
        print(f"Results saved in: {self.output_directory}")
        print(f"{'='*60}")


def main():
    """Main function to run the spectral analysis"""

    # Configuration
    catalog_file = "/orange/adamginsburg/salt/dihca2imaging/dihca_source_catalog.fits"
    image_directory = "/orange/adamginsburg/salt/dihca2imaging/grouped_imaging_products"
    output_directory = "/orange/adamginsburg/salt/dihca2imaging/spectral_analysis_results"

    # Initialize analyzer
    analyzer = SpectralAnalyzer(catalog_file, image_directory, output_directory)

    print("Starting spectral analysis of DIHCA sources...")

    return analyzer


def _get_field_from_filename(filename: str) -> str:
    """Extract the source_field prefix from a filename.

    Expected filename format: <source_field>_group..._spw... .image.fits
    """
    base = os.path.basename(filename)
    if '_group' in base:
        return base.split('_group', 1)[0]
    # fallback: take prefix before first underscore
    return base.split('_', 1)[0]


def _normalize_path(p: str) -> str:
    return os.path.abspath(os.path.expanduser(p))


def process_cube_cli(analyzer: SpectralAnalyzer, cube_path: str):
    """Process a single cube file: find catalog sources in the same field and run analysis only for that cube."""
    cube_path = _normalize_path(cube_path)
    if not os.path.exists(cube_path):
        raise SystemExit(f"Cube file not found: {cube_path}")

    field = _get_field_from_filename(cube_path)
    print(f"Processing cube {cube_path} for field {field}")

    # Find matching sources in the catalog
    matches = [row for row in analyzer.catalog if row['source_field'] == field]
    if not matches:
        print(f"No catalog sources found for field {field}; nothing to do.")
        return

    for row in matches:
        print(f"Running analysis for source {row['source_id']} using cube {os.path.basename(cube_path)}")
        analyzer.analyze_source(row, image_files=[cube_path])


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run spectral analysis or process a single cube as a job")
    parser.add_argument("--process-cube", help="Path to a single image cube to process (one cube per job)")
    parser.add_argument("--max-sources", type=int, default=None, help="Limit number of sources (when running full analysis)")
    args = parser.parse_args()

    analyzer = main()

    if args.process_cube:
        process_cube_cli(analyzer, args.process_cube)
    else:
        # For testing, analyze just the first few sources
        print("\nRunning test analysis on first 3 sources...")
        analyzer.analyze_all_sources(max_sources=args.max_sources or 3)
