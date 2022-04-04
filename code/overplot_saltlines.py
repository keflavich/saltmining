from spectral_cube import SpectralCube
import numpy as np
import regions
import os
import glob
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord
import radio_beam
import pyspeckit
import pylab as pl
from astroquery.splatalogue.utils import minimize_table as mt
import sys


sys.path.append('/orange/adamginsburg/salt/Orion_ALMA_2016.1.00165.S/analysis')
import lines
import paths

from salt_tables import (salt_tables, salt_table_names, SO, SO2, HCl, sis_tables, AlCl, AlF, Al37Cl,
                         NaF, AlO, AlOH, NaCN, CaS, CaO)


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

tables = salt_tables

def linename(row):
    return row['Species'] #"{0} {1}".format(row['Species'], row['QNs'])
def freq(row):
    return u.Quantity(row['Freq'], u.GHz)

linenames = [linename(row) for tbl in tables for row in tbl if row['vu'] <=4]
linetexnames = [linename(row) for tbl in tables for row in tbl if row['vu'] <=4] + ided_linetexnames
linetexnames = np.array(linetexnames)
linefreqs = np.hstack([u.Quantity([freq(row) for tbl in tables for row in tbl if row['vu'] <=4], u.GHz).value,
                       ided_linefreqs.value])
linefreqs = u.Quantity(linefreqs, u.GHz)


detection_table = table.Table.read(paths.tpath('salts_in_band.ipac'), format='ascii.ipac')
nondetections = (detection_table['Flag'] == '-n') | (detection_table['Flag'] == 'cn')
detection_table = detection_table[~nondetections]


def overplot_saltlines(spectra, vcen = 0*u.km/u.s, savepath='.', ymax=None):
    # (vcen in the stacked centroid is defined to be zero)

    fpath = lambda x: f"{savepath}/{x}"


    for sp_st in spectra:

        sp_st.xarr.convert_to_unit(u.GHz)

        basefn = os.path.basename(sp_st.specname)
        print(basefn)

        pl.figure(0, figsize=(16,6)).clf()
        sp_st.plotter(figure=pl.figure(0, figsize=(16,6)), clear=True)

        lines_to_plot = ((linefreqs > sp_st.xarr.as_unit(linefreqs.unit).min()*(1-vcen/constants.c)) &
                         (linefreqs < sp_st.xarr.as_unit(linefreqs.unit).max()*(1+vcen/constants.c)))
        
        if ymax is None:
            ymax = sp_st.data.max()
            if ymax < 0.004:
                ymax = 0.004
        sp_st.plotter(ymax=ymax)

        sp_st.plotter.line_ids(linetexnames[lines_to_plot], linefreqs[lines_to_plot], velocity_offset=vcen,
                               label1_size=16,
                               auto_yloc_fraction=0.75)
        for txt in sp_st.plotter.axis.texts:
            txt.set_backgroundcolor((1,1,1,0.9))


        sp_st.plotter.savefig(fpath('{0}'.format(basefn.replace("fits","png"))))

        for obj in sp_st.plotter.axis.texts+sp_st.plotter.axis.lines:
            if 'Na' in obj.get_label():
                obj.set_color('r')
                obj.set_zorder(5)
            elif 'K' in obj.get_label():
                obj.set_color('b')
                obj.set_zorder(10)
        for txt in sp_st.plotter.axis.texts:
            txt.set_backgroundcolor((1,1,1,0.9))

        sp_st.plotter.savefig(fpath('color_labels_{0}'
                                          .format(basefn.replace("fits","png")))
                              )

if __name__ == '__main__':
    basepath = '/orange/adamginsburg/salt/'
    from archive_metadata import field_meta

    import glob

    for field, field_data in field_meta.items():
        files = field_data['files']
        regions = regions.Regions.read(f'{basepath}/archive/regions/{field_data["regions"]}')
        for fn in files:
            #cube = SpectralCube.read(f'{basepath}/archive/{fn}')

            for reg in regions:
                name = reg.meta['text']
                basename = os.path.basename(fn)

                sp = pyspeckit.Spectrum(f'{basepath}/archive/{field}/spectra/Source{name}_{basename}')
                sp.specname = f'Source{name}_{basename}'

                if not os.path.exists(f'{basepath}/archive/{field}/spectra/figures'):
                    os.mkdir(f'{basepath}/archive/{field}/spectra/figures')

                overplot_saltlines([sp], vcen=u.Quantity(field_data['velocity'][name], u.km/u.s),
                                   savepath=f'{basepath}/archive/{field}/spectra/figures')
