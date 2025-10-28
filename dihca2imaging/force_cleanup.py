import shutil, os, glob, numpy as np
from spectral_cube import SpectralCube
from spectral_cube.utils import SpectralCubeWarning
import warnings
import subprocess


sacct = subprocess.check_output([
    '/opt/slurm/bin/sacct',
    '--format=JobID,JobName%60,Account%15,QOS%17,State,Priority%8,ReqMem%8,CPUTime%15,Elapsed%15,Timelimit%15,NodeList%20'
])

# Parse the output to get job information
lines = sacct.decode().split('\n')
jobs = {}
for line in lines[2:]:  # Skip header lines
    if line.strip():
        parts = line.split()
        if len(parts) >= 5:
            job_id = parts[0]
            job_name = parts[1]
            state = parts[4]
            jobs[job_name] = {'job_id': job_id, 'job_name': job_name, 'state': state}

def check_job_status(job_name):
    # use the directory name to avoid dealing with the .128 stuff
    #job_name = os.path.basename(os.path.split(image_name)[0])
    assert len(job_name) > 0
    for jn, job in jobs.items():
        if job_name in jn or jn in job_name or job['job_id'] in job_name or job_name in job['job_id']:
            return job
    return None

for fn in glob.glob("//red/adamginsburg/dihca/workdir_grouped/*/"):
    job_name = os.path.basename(fn.strip("/"))
    job_state = check_job_status(job_name)
    if job_state is not None and job_state['state'] in ('RUNNING', 'PENDING'):
        print(f"{fn} is still running, skipping ({job_state})")
        continue

    else:
        bad_files = glob.glob(os.path.join(fn, "IMAGING_WEIGHT*"))
        for bad_file in bad_files:
            print(f"removing {bad_file}")
            shutil.rmtree(bad_file)
        bad_files = glob.glob(os.path.join(fn, "TempLattice*"))
        for bad_file in bad_files:
            print(f"removing {bad_file}")
            shutil.rmtree(bad_file)


"""
Check the final products, then delete the /red workdir if the products are all there and check out
"""
for fn in glob.glob("/orange/adamginsburg/salt/dihca2imaging/grouped_imaging_products/*image.pbcor.fits"):

    job_name = os.path.basename(os.path.split(fn)[0])
    job_state = check_job_status(job_name)

    if job_state is not None and job_state['state'] in ('RUNNING', 'PENDING'):
        print(f"{fn} is still running, skipping ({job_state})")
        continue

    if not os.path.exists(fn):
        print(f"{fn} didn't exist")
        continue
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=SpectralCubeWarning)
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        cube = SpectralCube.read(fn, use_dask=True)
        mid = cube.shape[0]//2
        OK = False
        if np.any(np.isfinite(cube[mid,:,:])):
            print(f"{fn} is OK")
            OK = True
        else:
            print(f"{fn} has a bad middle channel")
            mx = cube.max()
            if np.isnan(mx):
                print(f"{fn} max={mx}, deleting")
                for suffix in ('image', 'model', 'pb', 'psf', 'residual', 'sumwt', 'wt', 'image.pbcor', 'mask'):
                    torm = fn.replace("image", suffix)
                    if os.path.exists(torm):
                        print(f"removing {torm}")
                        shutil.rmtree(torm)
                    if os.path.exists(torm+".fits"):
                        print(f"removing {torm}.fits")
                        os.remove(torm+".fits")
            else:
                OK = True
                print(f"{fn} is OK")

    if OK:
        base = os.path.basename(fn).replace(".image.pbcor.fits", "").replace(".image.fits", "").replace("image.pbcor","").replace(".image", "")
        torm = os.path.join('/red/adamginsburg/dihca/workdir_grouped/', base)
        if os.path.exists(torm):
            print(f"removing {torm}")
            shutil.rmtree(torm)
        else:
            print(f"{torm} doesn't exist - we're totally done with it!")

def has_beam(cube):
    try:
        cube.beam
    except AttributeError:
        try:
            cube.beams
        except AttributeError:
            return False
    return True

"""
Check the individual images
"""
for fn in glob.glob("/red/adamginsburg/dihca/workdir_grouped/*/*image"):
    job_name = os.path.basename(os.path.split(fn)[0])
    job_state = check_job_status(job_name)
    if job_state is not None and job_state['state'] in ('RUNNING', 'PENDING'):
        print(f"{fn} is still running, skipping ({job_state})")
        continue

    if not os.path.exists(fn):
        print(f"{fn} didn't exist")
        continue
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=SpectralCubeWarning)
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        cube = SpectralCube.read(fn, use_dask=True)
    mid = cube.shape[0]//2
    isfin = np.any(np.isfinite(cube[mid,:,:]))
    if isfin and has_beam(cube):
        print(f"{fn} is OK")
        continue
    else:
        if not has_beam(cube):
            print(f"{fn} has no beam")
        if not isfin:
            print(f"{fn} has a bad middle channel")
        mx = cube.max()
        if np.isnan(mx):
            print(f"{fn} max={mx}, deleting")
            for suffix in ('image', 'model', 'pb', 'psf', 'residual', 'sumwt', 'wt', 'image.pbcor', 'mask'):
                torm = fn.replace("image", suffix)
                if os.path.exists(torm):
                    print(f"removing {torm}")
                    shutil.rmtree(torm)
                if os.path.exists(torm+".fits"):
                    print(f"removing {torm}.fits")
                    os.remove(torm+".fits")
        else:
            print(f"{fn} is OK")