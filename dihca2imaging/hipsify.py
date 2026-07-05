import numpy as np
import os
import shutil
import PIL.Image
import pyavm
from astropy import log
from scipy.ndimage import label, binary_dilation

from tqdm.auto import tqdm
from reproject.hips import reproject_to_hips
from reproject import reproject_interp

import glob

log.setLevel('INFO')

def convert_black_to_transparent(image_path):
    """Convert solid black pixels to transparent in an image, but only if they're connected to the edges"""
    img = PIL.Image.open(image_path)

    # Convert to RGBA if not already
    if img.mode != 'RGBA':
        img = img.convert('RGBA')

    # Convert to numpy array for easier manipulation
    img_array = np.array(img)
    height, width = img_array.shape[:2]

    # Find solid black pixels (R=0, G=0, B=0)
    black_pixels = (img_array[:, :, 0] == 0) & (img_array[:, :, 1] == 0) & (img_array[:, :, 2] == 0)

    # Label connected components of black pixels
    labeled_regions, num_regions = label(black_pixels)

    # Create a mask for edge-connected black pixels
    edge_connected_black = np.zeros_like(black_pixels, dtype=bool)

    # Check which labeled regions touch the edges
    edge_touching_labels = set()

    # Check top and bottom edges
    edge_touching_labels.update(labeled_regions[0, :])  # Top edge
    edge_touching_labels.update(labeled_regions[height-1, :])  # Bottom edge

    # Check left and right edges
    edge_touching_labels.update(labeled_regions[:, 0])  # Left edge
    edge_touching_labels.update(labeled_regions[:, width-1])  # Right edge

    # Remove label 0 (background/non-black pixels)
    edge_touching_labels.discard(0)

    # Create mask for all pixels belonging to edge-touching regions
    for label_id in edge_touching_labels:
        edge_connected_black |= (labeled_regions == label_id)

    # Set alpha to 0 only for edge-connected black pixels
    img_array[edge_connected_black, 3] = 0

    # Convert back to PIL Image
    img_transparent = PIL.Image.fromarray(img_array, 'RGBA')

    # Save with transparent suffix
    transparent_path = image_path.replace('.png', '_transparent.png').replace('.jpg', '_transparent.jpg')
    if transparent_path.endswith('_transparent.jpg'):
        transparent_path = transparent_path.replace('_transparent.jpg', '_transparent.png')
    img_transparent.save(transparent_path)

    avm = pyavm.AVM.from_image(image_path)
    avm.embed(transparent_path, transparent_path)

    return transparent_path

def fits_to_png(filename):
    from astropy.visualization import simple_norm
    import matplotlib.colors as mcolors
    from astropy.io import fits
    import pylab as pl
    from astropy.wcs import WCS
    import pyavm
    import PIL

    colors1 = pl.cm.gray_r(np.linspace(0., 1, 128))
    colors2 = pl.cm.hot(np.linspace(0, 1, 128))

    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    data = fits.getdata(filename).squeeze()
    #normed = simple_norm(data, stretch='asinh', max_percent=99.995, min_percent=1)(data)
    #colorized = (mymap(normed) * 255).astype('uint8')

    # try a bilinear approach
    std = np.nanstd(data)
    normed1 = simple_norm(data, stretch='linear', vmin=-5*std, vmax=5*std)(data)
    normed2 = simple_norm(data, stretch='sqrt', vmin=5*std)(data) # vmax=max
    colorized1 = (pl.cm.gray_r(normed1) * 255).astype('uint8')
    colorized2 = (pl.cm.hot(normed2) * 255).astype('uint8')
    colorized = colorized1
    mask = data > 5*std
    colorized[mask, :] = colorized2[mask, :]

    # flip along y-axis
    img = PIL.Image.fromarray(colorized[::-1, :, :])
    pngfn = filename.replace('fits', 'png')
    img.save(pngfn)

    wcs = WCS(fits.getheader(filename)).celestial

    AVM = pyavm.AVM.from_wcs(wcs, shape=data.shape)
    AVM.embed(pngfn, pngfn)

    return pngfn


def main():
    basepath = '/orange/adamginsburg/salt/dihca/'

    import PIL
    import PIL.Image
    PIL.Image.MAX_IMAGE_PIXELS = 207360000 * 4

    for globstr, savename in [
                              (f'{basepath}/*fits', 'DIHCA2_continuum_HIPS'),
                              ]:
        hips_list = []
        filelist = glob.glob(globstr)
        print(f"Processing {len(filelist)} files for {savename}:")
        print(filelist, flush=True)
        for filename in filelist:
            print(filename, flush=True)
            if filename.endswith('fits'):
                #filename = fits_to_png(filename)
                # switch to this if you _don't_ want to remake the pngs
                if not os.path.exists(filename.replace("fits", "png")):
                    filename = fits_to_png(filename)
                else:
                    filename = filename.replace("fits", "png")
            try:
                avm = pyavm.AVM.from_image(filename)
            except pyavm.exceptions.NoXMPPacketFound:
                print(f"No XMP packet found for {filename}")
                continue
            except OSError:
                os.remove(filename)
                filename = fits_to_png(filename.replace("png", "fits"))
                avm = pyavm.AVM.from_image(filename)

            # Convert black pixels to transparent
            # print(f"Converting black pixels to transparent for {filename}")
            # filename_transparent = filename.replace('.png', '_transparent.png').replace('.jpg', '_transparent.jpg')
            # if not os.path.exists(filename_transparent):
            #     filename_transparent = convert_black_to_transparent(filename)

            # no need to make black transparent
            filename_transparent = filename

            # Copy AVM metadata to the transparent version
            avm.embed(filename_transparent, filename_transparent)

            # Use the transparent version for processing
            processing_filename = filename_transparent

            output_directory = processing_filename.replace('.png', '_hips').replace('.jpg', '_hips')

            if os.path.exists(output_directory):
                # Compare modification times: if the source image is newer than the
                # existing hips directory (or any file inside it), remove the
                # directory so it will be remade. If the directory is newer or
                # equal, skip and add to the list.
                file_mtime = os.path.getmtime(processing_filename)

                # Find the newest modification time inside the output directory
                latest_dir_mtime = os.path.getmtime(output_directory)
                for root, dirs, files in os.walk(output_directory):
                    for fname in files:
                        try:
                            p = os.path.join(root, fname)
                            mt = os.path.getmtime(p)
                            if mt > latest_dir_mtime:
                                latest_dir_mtime = mt
                        except OSError:
                            # ignore unreadable files
                            continue

                if file_mtime <= latest_dir_mtime:
                    # The existing hips folder is up-to-date (or newer) than the
                    # source image; skip reprocessing.
                    print(f"Found up-to-date directory; skipping {filename}")
                    hips_list.append(output_directory)
                    continue
                else:
                    # Source image is newer than the existing hips output; remove
                    # the folder so it will be remade by reproject_to_hips.
                    print(f"Source {processing_filename} is newer than {output_directory}; removing and remaking.")
                    shutil.rmtree(output_directory)

            print("filename, processing_filename, output_directory: ", filename, processing_filename, output_directory, )#np.array(PIL.Image.open(processing_filename)).shape ) #pyavm.AVM.from_image(filename))


            reproject_to_hips(processing_filename,
                        coord_system_out='galactic',
                        level=None,
                        reproject_function=reproject_interp,
                        output_directory=output_directory,
                        threads=16,
                        progress_bar=tqdm)

            hips_list.append(output_directory)

        from reproject.hips import coadd_hips
        if os.path.exists(f'{basepath}/{savename}'):
            shutil.rmtree(f'{basepath}/{savename}')
        assert len(hips_list) > 0
        coadd_hips(hips_list, f'{basepath}/{savename}')

if __name__ == "__main__":
    main()