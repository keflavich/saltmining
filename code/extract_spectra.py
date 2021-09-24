import regions
from spectral_cube import SpectralCube
import os
from archive_metadata import field_meta
basepath = '/orange/adamginsburg/salt/'

for field, field_data in field_meta.items():
    files = field_data['files']
    regions = regions.Regions.read(f'{basepath}/archive/regions/{field_data["regions"]}')
    for fn in files:
        cube = SpectralCube.read(f'{basepath}/archive/{fn}')

        for reg in regions:
            subcube = cube.subcube_from_regions([reg])
            spec = subcube.mean(axis=(1,2))
            name = reg.meta['text']
            basename = os.path.basename(fn)

            if not os.path.exists(f'{basepath}/archive/{field}/'):
                os.mkdir(f'{basepath}/archive/{field}/')
            if not os.path.exists(f'{basepath}/archive/{field}/spectra/'):
                os.mkdir(f'{basepath}/archive/{field}/spectra/')

            spec.write(f'{basepath}/archive/{field}/spectra/Source{name}_{basename}')
