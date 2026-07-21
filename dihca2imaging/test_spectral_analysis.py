# test_spectral_analysis.py

import time
from spectral_cube import SpectralCube
import sys
sys.path.append('/orange/adamginsburg/salt/dihca2imaging/')
import spectral_analysis
from astropy.coordinates import SkyCoord
from astropy import units as u

analyzer = spectral_analysis.main()

analyzer.catalog.add_index('source_field')
analyzer.catalog.add_index('source_number')

row = analyzer.catalog.loc['G012.908-00.259'].loc[('source_number', 5)]
source_coord = SkyCoord(row['ra_deg'], row['dec_deg'], unit=(u.deg, u.deg))

cube = SpectralCube.read('/orange/adamginsburg/salt/dihca2imaging/grouped_imaging_products/G012.908-00.259_group2_spw0.image.fits', use_dask=True)
cube.allow_huge_operations = True
print("Computing mad_std of first 50 channels: ", end="", flush=True)
t0 = time.time()
noise_level = cube[:50].mad_std()
print(f"{time.time()-t0:0.3f} seconds")
mask = analyzer.create_source_mask(cube, source_coord, noise_level)

spectrum = analyzer.extract_spectrum(cube, mask)

# Step 3: Find spectral peaks
peaks_catalog = analyzer.find_spectral_peaks(spectrum, noise_level)

consolidated_peaks = analyzer.consolidate_peaks(peaks_catalog, velocity_threshold=5.0*u.km/u.s)

# Step 4: Determine source velocity (using consolidated peaks)
velocity, velocity_method = analyzer.determine_source_velocity(source_coord, consolidated_peaks)

# Step 4.5: Attempt to identify the remaining line peaks
consolidated_peaks = analyzer.identify_line_peaks(consolidated_peaks, velocity)

# Step 5: Create moment images for each consolidated peak
all_moments = analyzer.create_moment_images_for_peaks(cube, mask, consolidated_peaks)

# Step 6: Create masked moment images for each consolidated peak
all_masked_moments = analyzer.create_masked_moment_images_for_peaks(cube, mask, consolidated_peaks)
