"""
Retrieve data from the ALMA archive and look for salt

The archive retrieval is only half-completed; some hacking is needed.

The archival image cubes only cover a small fraction of the bandwidth observed.
"""
import numpy as np
import scipy.ndimage
from astroquery.alma import Alma
from astropy import units as u, constants
from spectral_cube import SpectralCube
import spectral_cube.analysis_utilities
import pyspeckit
from astropy import stats
#import reproject
import os

import pylab as pl

if not os.path.exists('member.uid___A001_X12a_X19.A2591_line_spw1.image.pbcor.fits'):
    rslt = Alma.query_object('AFGL2591')
    data = Alma.retrieve_data_from_uid(rslt['member_ous_id'])


cube = SpectralCube.read('member.uid___A001_X12a_X19.A2591_line_spw1.image.pbcor.fits')

# many from https://www.aanda.org/articles/aa/pdf/2019/11/aa35865-19.pdf
# but the ?'s were *not detected* in that article!  What the hell? (beam was 1.5")
known_lines = {'13CS': 231.220996*u.GHz,
               'CH3OH1029': 231.28115*u.GHz,
               'CH3CO?': 231.3104984*u.GHz,
               'H213CO': 231.2459639*u.GHz,
               'HNO3b?': 231.2517153*u.GHz,
               '(CH3)2COd?': 231.1999571*u.GHz,
               '34SO2': 219.3550091*u.GHz,
               'SO2v2=1': 219.466*u.GHz,
               'SO2': 219.276*u.GHz,
               'HC3Nv=1': 219.174*u.GHz,
               'HNO3?': 219.3838408*u.GHz,
               'MnO??': 219.3045193*u.GHz,
               '(CH3)2COa?': 219.2421408*u.GHz,
               '(CH3)2COb?': 219.2642749*u.GHz,
               '(CH3)2COc?': 219.219931*u.GHz,
              }




mx = cube.spectral_slab(231.210*u.GHz, 231.220*u.GHz)[:, 460:560, 480:580].max(axis=0)
mask = scipy.ndimage.binary_dilation((mx > stats.mad_std(mx)*7), iterations=6)

vmap = cube.spectral_slab(231.210*u.GHz, 231.220*u.GHz)[:, 460:560, 480:580].with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=231.220996*u.GHz).argmax_world(axis=0)

vmap[~mask] = np.nan
#vmap_masked = np.ma.masked_where(~mask, vmap.value)


#for fullcube in (cube[:, 460:560, 480:580], ):
for spw in (1,2):
#for fn in ('member.uid___A001_X12a_X19.A2591_line_spw1.image.pbcor.fits', 'member.uid___A001_X12a_X19.A2591_line_spw2.image.pbcor.fits'):
    fn = f'member.uid___A001_X12a_X19.A2591_line_spw{spw}.image.pbcor.fits'
    fullcube = SpectralCube.read(fn)[:, 460:560, 480:580]
    # convert the cube to velocity units with an arbitrary reference point
    # (this step assumes the cube is in frequency or wavelength; if the
    # cube is not, it should be skipped)
    fullcube = fullcube.with_spectral_unit(u.km/u.s,
                                           velocity_convention='radio',
                                           rest_value=fullcube.spectral_axis.mean())

    # mask out super bright SiO masers; they break the FFT shifting tool
    # (this step can be skipped if there's nothing anomalously bright
    # in your spectrum)
    # fullcube = fullcube.with_mask(fullcube < 0.5*u.Jy/u.beam)

    # reproject the velocity map into the cube's coordinate system
    # * NOT NEEDED*
    #vmap_proj,_ = reproject.reproject_interp(vmap.hdu,
    #                                         fullcube.wcs.celestial,
    #                                         shape_out=fullcube.shape[1:])
    #vmap_proj = u.Quantity(vmap_proj, u.km/u.s)
    vmap_proj = vmap

    # perform the stacking!
    stack = spectral_cube.analysis_utilities.stack_spectra(fullcube, vmap_proj,
                                                           v0=0.0*u.km/u.s)
    fstack = stack.with_spectral_unit(u.GHz)

    fstack.write(f'stacked_spectra/AFGL2591_stack_spw{spw}.fits', overwrite=True)

    pl.clf()
    fstack.quicklook(filename=f'stacked_spectra/AFGL2591_stack_spw{spw}.PDF')


import numpy as np
import os
import glob
import paths
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord
import radio_beam
import pyspeckit
import pylab as pl
from astroquery.splatalogue import Splatalogue
from astroquery.splatalogue.utils import minimize_table as mt
import sys
sys.path.append('../../orion/analysis')
import lines
from salt_tables import (salt_tables, salt_table_names, SO, SO2, HCl, sis_tables, AlCl, AlF, Al37Cl,
                         NaF, AlO, AlOH, NaCN, CaS, CaO)

# TO DETERMINE:
# (vcen in the stacked centroid is defined to be zero)
vcen = 0.*u.km/u.s

pl.matplotlib.rcParams['font.size']=16

all_lines = {**lines.disk_lines, **lines.absorbers}

ided_linenames = sorted(all_lines.keys())
ided_linefreqs = u.Quantity([all_lines[x] for x in ided_linenames
                             #if 'U' not in x
                            ])
ided_linetexnames = [lines.texnames[x] if x in lines.texnames else x
                     for x in ided_linenames
                     #if 'U' not in x
                    ]

#salt_tables = [KCl, K37Cl, K41Cl, NaCl, Na37Cl, K41Cl37]
salt_colors = ['b', 'm', 'darkgreen', 'orange', 'c', 'y']

tables = []

def linename(row):
    return "{0} {1}".format(row['Species'], row['QNs'])
def freq(row):
    return u.Quantity(row['Freq'], u.GHz)

linenames = [linename(row) for tbl in tables for row in tbl]
linetexnames = [linename(row) for tbl in tables for row in tbl] + ided_linetexnames
linefreqs = np.hstack([u.Quantity([freq(row) for tbl in tables for row in tbl], u.GHz).value,
                       ided_linefreqs.value])
linefreqs = u.Quantity(linefreqs, u.GHz)

detection_table = table.Table.read(paths.tpath('salts_in_band.ipac'), format='ascii.ipac')
nondetections = (detection_table['Flag'] == '-n') | (detection_table['Flag'] == 'cn')
detection_table = detection_table[~nondetections]


paths.fpath = lambda x: x


flist = glob.glob('stacked_spectra/AFGL2591*.fits')
for fn in flist:
    print(fn)

    sp_st = pyspeckit.Spectrum(fn)

    pl.figure(0, figsize=(16,6)).clf()
    sp_st.plotter(figure=pl.figure(0, figsize=(16,6)), clear=True)

    basefn = os.path.split(fn)[-1]

    sp_st.plotter(ymin=-0.0025, ymax=0.01)
    sp_st.plotter.line_ids(linetexnames, linefreqs, velocity_offset=vcen,
                           label1_size=16,
                           auto_yloc_fraction=0.75)
    for txt in sp_st.plotter.axis.texts:
        txt.set_backgroundcolor((1,1,1,0.9))


    sp_st.plotter.savefig(paths.fpath('stacked_spectra/{0}'
                                      .format(basefn.replace("fits","pdf")))
                          )

    for obj in sp_st.plotter.axis.texts+sp_st.plotter.axis.lines:
        if 'Na' in obj.get_label():
            obj.set_color('r')
            obj.set_zorder(5)
        elif 'K' in obj.get_label():
            obj.set_color('b')
            obj.set_zorder(10)
    for txt in sp_st.plotter.axis.texts:
        txt.set_backgroundcolor((1,1,1,0.9))

    sp_st.plotter.savefig(paths.fpath('stacked_spectra/color_labels_{0}'
                                      .format(basefn.replace("fits","pdf")))
                          )


    sp_st.plotter(ymin=-0.0025, ymax=0.01)

    # uses lines.py
    sp_st.plotter.line_ids(linetexnames, linefreqs, velocity_offset=vcen,
                           auto_yloc_fraction=0.8)

    for txt in sp_st.plotter.axis.texts:
        txt.set_backgroundcolor((1,1,1,0.9))
    #sp_st.plotter.line_ids(ided_linetexnames, ided_linefreqs, velocity_offset=-vcen,
    #                       plot_kwargs=dict(color='b'))
    sp_st.plotter.savefig(paths.fpath('stacked_spectra/lines_labeled_{0}'
                                      .format(basefn.replace("fits","pdf")))
                         )

    #sp_st.plotter(ymin=-0.0025, ymax=0.01)
    # use the salt names directly.  This is for labeling of the colored
    # lines; the publication-ready stuff still uses lines.py
    #sp_st.plotter.line_ids(detection_table['Species'],
    #                       u.Quantity(detection_table['Frequency'], u.GHz),
    #                       velocity_offset=-vcen,
    #                       auto_yloc_fraction=0.8)

    for tbl,color,nm in zip(salt_tables, salt_colors, salt_table_names):
        for row in tbl:
            frq = u.Quantity(row['Freq'], u.GHz)
            if frq > sp_st.xarr.min() and frq < sp_st.xarr.max():
                print(row)
                sp_st.plotter.axis.vlines((frq*(1-vcen/constants.c)).to(u.GHz).value,
                                          -0.05, 0.10,
                                          colors=color, linestyles=':', label=nm+row['QNs'])

    for linename, linefreq in known_lines.items():
        sp_st.plotter.axis.vlines((linefreq*(1-vcen/constants.c)).to(u.GHz).value,
                                  -0.05, 0.10,
                                  colors='g' if '?' in linename else 'r', linestyles='--', label=linename)

    #for row in HCl:
    #    frq = u.Quantity(row['Freq'], u.GHz).value
    #    if frq > sp_st.xarr.min().value and frq < sp_st.xarr.max().value:
    #        sp_st.plotter.axis.vlines(frq*(1-vcen/constants.c).decompose().value,
    #                                  -0.05, 0.10,
    #                                  colors='g', linestyles='--')
    pl.legend(loc='upper right')

    sp_st.plotter.savefig(paths.fpath('stacked_spectra/diagnostic_lines_labeled_{0}'
                                      .format(basefn.replace("fits","pdf")))
                         )

    # Do another one just for SiO
    sp_st.plotter(ymin=-0.0025, ymax=0.01)

    # uses lines.py
    sp_st.plotter.line_ids(linetexnames, linefreqs, velocity_offset=vcen,
                           auto_yloc_fraction=0.8)

    for txt in sp_st.plotter.axis.texts:
        txt.set_backgroundcolor((1,1,1,0.9))


    for tbl,color in zip(sis_tables, salt_colors):
        for row in tbl:
            frq = u.Quantity(row['Freq'], u.GHz).value
            if frq > sp_st.xarr.min().value and frq < sp_st.xarr.max().value:
                sp_st.plotter.axis.vlines(frq*(1-vcen/constants.c).decompose().value,
                                          -0.05, 0.10,
                                          colors=color, linestyles=':')

    sp_st.plotter.savefig(paths.fpath('stacked_spectra/sis_diagnostic_lines_labeled_{0}'
                                      .format(basefn.replace("fits","pdf")))
                         )

    # Do another one just for alcl
    sp_st.plotter(ymin=-0.0025, ymax=0.01)

    # uses lines.py
    sp_st.plotter.line_ids(linetexnames, linefreqs, velocity_offset=vcen,
                           auto_yloc_fraction=0.8)

    for txt in sp_st.plotter.axis.texts:
        txt.set_backgroundcolor((1,1,1,0.9))

    for tbl,color in zip([AlCl, AlF, Al37Cl, NaF, AlOH, AlO], salt_colors):
        for row in tbl:
            frq = u.Quantity(row['Freq'], u.GHz).value
            if frq > sp_st.xarr.min().value and frq < sp_st.xarr.max().value:
                sp_st.plotter.axis.vlines(frq*(1-vcen/constants.c).decompose().value,
                                          -0.05, 0.10,
                                          colors=color, linestyles=':')


    sp_st.plotter.savefig(paths.fpath('stacked_spectra/alcl_diagnostic_lines_labeled_{0}'
                                      .format(basefn.replace("fits","pdf")))
                         )


    for (a,b) in zip(linetexnames, linefreqs):
        if (b>sp_st.xarr.min()) and (b<sp_st.xarr.max()) and a not in ided_linetexnames:
            print("'{0}': {1}*u.{2},".format(a,b.value,b.unit))

    for speciesname, species in (('NaCN', NaCN), ('SO2',SO2), ('SO', SO), ('CaS', CaS), ('CaO', CaO)):
        # Do another one just for nacn
        sp_st.plotter(ymin=-0.0025, ymax=0.01)

        # uses lines.py
        sp_st.plotter.line_ids(linetexnames, linefreqs, velocity_offset=vcen,
                               auto_yloc_fraction=0.8)

        for txt in sp_st.plotter.axis.texts:
            txt.set_backgroundcolor((1,1,1,0.9))

        for tbl,color in zip([species], ['b']):
            for row in tbl:
                frq = u.Quantity(row['Freq'], u.GHz).value
                if frq > sp_st.xarr.min().value and frq < sp_st.xarr.max().value:
                    sp_st.plotter.axis.vlines(frq*(1-vcen/constants.c).decompose().value,
                                              -0.05, 0.10,
                                              colors=color, linestyles=':')


        sp_st.plotter.savefig(paths.fpath('stacked_spectra/{1}_diagnostic_lines_labeled_{0}'
                                          .format(basefn.replace("fits","pdf"),
                                                  speciesname
                                                 ))
                             )
