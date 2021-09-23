import os
from find_salt_obs import result, freqs
from astropy import units as u
from astroquery.alma import Alma, utils as almautils
from astroquery.splatalogue import Splatalogue, utils as splatutils
from spectral_cube import SpectralCube

if not os.path.exists('candidates'):
    os.mkdir('candidates')

for row in result:
    if row['Has salt']:

        sourcename = row['Source name']

        mous = row['Member ous id']
        reg_mous = mous.replace(":","_").replace("/","_")
        if os.path.exists(reg_mous):
            continue

        stage = Alma().stage_data(mous)

        tarfiles = [x for x in stage['URL'] if '001_of_001.tar' in x]
        fitsfiles = [x for x in stage['URL'] if 'cube.I.pbcor.fits' in x]

        almadl = Alma()
        almadl.cache_location = '.'

        if len(fitsfiles) > 0:
            localfiles = almadl.download_and_extract_files(fitsfiles)
        else:
            print("SKIP")
            print(fitsfiles, tarfiles)
            continue
            localfiles = almadl.download_and_extract_files(tarfiles)

        for fn in localfiles:
            print(fn)
            cube = SpectralCube.read(fn)
            lines = cube.find_lines(chemical_name='NaCl|KCl', line_lists=['CDMS'])
            if any(lines):
                lines = splatutils.minimize_table(lines)

            for line in lines:
                print(line)
                subcube = (cube
                           .with_spectral_unit(u.km/u.s,
                                               velocity_convention='radio',
                                               rest_value=u.Quantity(line['Freq'],
                                                                     u.GHz))
                           .spectral_slab(-100*u.km/u.s, 100*u.km/u.s)
                          )
                if len(subcube) > 1:
                    linename = "{0}_{1}_{2}".format(line['Species'],
                                                    line['QNs'],
                                                    line['Freq'])
                    subcube.write(f'candidates/{sourcename}_{linename}_cube.fits', overwrite=True)
                    subcube.max(axis=0, how='slice').write(f'candidates/{sourcename}_{linename}_max.fits', overwrite=True)

        with open(reg_mous, 'w') as fh:
            fh.write("DONE")
