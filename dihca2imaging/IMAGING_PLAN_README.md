# DIHCA2 Imaging Plan

This document describes the automatically generated CASA imaging plan for the DIHCA2 source catalog.

## Overview

The imaging plan groups 533 sources from 65 fields into 211 imaging groups, resulting in 844 total imaging commands (211 groups × 4 spectral windows).

## Methodology

### Source Grouping

Sources were grouped using hierarchical clustering with the following criteria:

1. **Spatial Clustering**: Sources within 5.12" of each other are grouped together (matching the default 512×512 pixel field of view at 0.01"/pixel)

2. **Adaptive Image Sizes**:
   - Default: 512×512 pixels (5.12" FOV) for small, compact groups
   - Medium: 1024×1024 pixels (10.24" FOV) for groups with 5+ sources spread over larger areas
   - Large: 1536×1536 pixels (15.36" FOV) for very extended groups (not needed in current catalog)

3. **Group Optimization**: The algorithm minimizes the number of images while ensuring all sources are covered with appropriate spatial resolution

### Imaging Parameters

Each imaging command includes:
- **vis**: Measurement set path (`/orange/adamginsburg/ALMA_IMF/Uploads/dihca2/original_uvdata/`)
- **field**: CASA field name (e.g., "008.672-00.682")
- **spw**: Spectral window ID (0-3)
- **phasecenter**: J2000 coordinates centered on the source group
- **imsize**: Image size in pixels [width, height]
- **cell**: Pixel size (0.01 arcsec)
- **nchan**: Number of channels (3840)
- **gridder**: "mosaic" for multi-pointing observations
- **specmode**: "cube" for spectral line imaging
- **deconvolver**: "hogbom" algorithm
- **weighting**: Briggs with robust=0.5
- **niter**: 10000 iterations
- **threshold**: 0.1 mJy stopping threshold

### Spectral Windows

Four spectral windows are imaged for each group:
- **SPW 0**: H₂CO 3(0,3)-2(0,2) at 218.541687 GHz
- **SPW 1**: ¹³CO 2-1 at 220.398684 GHz
- **SPW 2**: SiO 5-4 at 217.104984 GHz
- **SPW 3**: C¹⁸O 2-1 at 219.560358 GHz

## Results Summary

### Statistics

- **Total imaging commands**: 844
- **Total groups**: 211
- **Total fields**: 65
- **Average groups per field**: 3.2

### Image Size Distribution

- **512×512 pixels (5.12" FOV)**: 159 groups (75%)
- **1024×1024 pixels (10.24" FOV)**: 52 groups (25%)
- **1536×1536 pixels (15.36" FOV)**: 0 groups (0%)

### Sources per Group

- **Minimum**: 1 source
- **Maximum**: 48 sources (G345.493+01.469_group3)
- **Mean**: 2.5 sources
- **Median**: 2 sources

### Top Fields by Number of Groups

1. **351.444+00.659**: 11 groups (40 sources)
2. **333.018+00.766**: 9 groups (26 sources)
3. **034.401+00.226**: 8 groups (35 sources)
4. **034.258+00.166**: 6 groups (19 sources)
5. **035.522-00.274**: 6 groups (7 sources)
6. **345.493+01.469**: 6 groups (54 sources)

## Output Files

### dihca2_imaging_plan.json

JSON file containing all imaging commands. Each entry has a unique key in the format:
```
{field_name}_group{N}_spw{M}
```

For example:
- `G008.672-00.682_group1_spw0`
- `G345.493+01.469_group3_spw2`

### Metadata

Each imaging command includes metadata:
- `n_sources`: Number of sources in this group
- `fov_arcsec`: Field of view in arcseconds
- `source_ids`: List of source IDs from the catalog

## Usage

To use this imaging plan in CASA:

```python
import json

# Load the imaging plan
with open('dihca2_imaging_plan.json', 'r') as f:
    imaging_plan = json.load(f)

# Execute a specific imaging command
cmd_key = 'G008.672-00.682_group1_spw0'
params = imaging_plan[cmd_key]

# Run tclean with the parameters
tclean(**params)
```

Or to process all commands:

```python
for cmd_key, params in imaging_plan.items():
    print(f"Imaging {cmd_key}...")
    # Remove metadata before passing to tclean
    tclean_params = {k: v for k, v in params.items() if k != 'metadata'}
    tclean(**tclean_params)
```

## Notes

- The measurement set paths assume data location at `/orange/adamginsburg/ALMA_IMF/Uploads/dihca2/original_uvdata/`
- For fields with `_eb1` suffix, the script attempts to find the base field MS file
- All images use the standard gridder appropriate for single-pointing ALMA observations
- Phase centers are computed as the mean position of all sources in each group
- A 20% padding is added around source extents to ensure full coverage

## Regenerating the Plan

To regenerate the imaging plan with different parameters, edit `imaging_plan.py` and run:

```bash
python imaging_plan.py
```

Key parameters to adjust:
- `max_separation_arcsec` in `cluster_sources()`: Controls how close sources must be to be grouped
- `max_size_arcsec` in `compute_group_bounds()`: Maximum allowed FOV (currently 15")
- `cell_size`: Pixel size in arcseconds (currently 0.01")
- Rest frequencies in `generate_tclean_commands()`: Spectral line frequencies for each SPW
