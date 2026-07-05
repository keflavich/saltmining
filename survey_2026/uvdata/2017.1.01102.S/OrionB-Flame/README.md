# OrionB-Flame / NGC 2024 Mosaic — 2017.1.01102.S

Downloaded 2026-06-29 via `download_orionb_flame.py`.

## What was delivered
- 1 sci continuum (mfs) FITS: `NGC_2024_Mosaic_sci.spw19_23_25_27_29.mfs.I.manual.image.pbcor.fits`
- 10 calibrator pbcor / pb fits (J0423-0120 bandpass, J0541-0211 phase)

## What's missing
- Spectral cube products (`*sci*cube*pbcor.fits`). The PI only delivered
  continuum imaging — to run the salt-search pipeline on this data we
  need to re-image from raw visibilities.

## Coverage
- Band 6: SPW 19, 23, 25, 27, 29 spanning 218.47, 219.54, 220.39, 230.59, 231.98 GHz
- These cover NaCl_v0/v1/v2 J=17-16 (217.97), J=18-17 (218.50), KCl_v0_J17-16 (217.78),
  H2O_5_15-4_22_232 (232.69)
- Beam: 0.142" → 58 AU at d=410 pc
- Useful for the salt search if cubes can be generated

## TODO to bring it into the pipeline
1. Request ASDM raw visibilities from ALMA archive (member_ous_uid = uid://A001/X129e/X541)
2. Run CASA tclean per-spw at native channel width with appropriate masking
3. Save under uvdata/2017.1.01102.S/OrionB-Flame/ as
   `member.<mous>.NGC_2024_Mosaic_sci.spw<NN>.cube.I.pbcor.fits`
4. Re-run launch_l4_d2_pipeline.py --force (will auto-pick this proposal)
