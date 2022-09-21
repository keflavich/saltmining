
def make_and_examine_moment_map(cubes, reg, vcen,
        pix_threshold=20*u.K,
        mxthreshold=20*u.K,
        erode_iter=1, dilate_iter=1,
        vrange=[-20,20]*u.km/u.s,
        restval=217.81764400*u.GHz):
    cube = [c for c in cubes if (c.spectral_extrema[0] < restval) & (c.spectral_extrema[1] > restval)][0]
    print(cube)
    scube = (cube.subcube_from_regions(reg)
             .with_spectral_unit(u.km/u.s,
                                 velocity_convention='radio',
                                 rest_value=restval)
             .spectral_slab(vrange[0]+vcen, vrange[1]+vcen))
    print(scube)
    pl.figure(figsize=(14,5))
    pl.subplot(1,3,1)
    m0 = scube.moment0()
    m0.quicklook()
    pl.title("m0")

    pl.subplot(1,3,2)
    mx = scube.max(axis=0)
    mx.quicklook()
    pl.title("mx")
    pl.colorbar()
    pl.subplot(1,3,3)
    msk = mx > mxthreshold
    if erode_iter > 0:
        msk = scipy.ndimage.binary_erosion(msk, iterations=erode_iter)
    if dilate_iter > 0:
        msk = scipy.ndimage.binary_dilation(msk, iterations=dilate_iter)
    pl.imshow(msk, interpolation='none', origin='lower')
    pl.title('mask')
    pl.tight_layout()
    pl.figure()


    m1 = scube.with_mask(msk).with_mask(scube > pix_threshold).moment1()
    pl.figure(figsize=(14,6))
    ax1 = pl.subplot(1,2,1)
    im = ax1.imshow(m1.value, vmin=-8, vmax=2)
    pl.colorbar(mappable=im)
    mxv = scube.with_mask(msk).argmax_world(axis=0)
    pl.subplot(1,2,2)
    pl.imshow(mxv.value, vmin=-8, vmax=2)
    pl.colorbar()

    vmap = m1.hdu

    return vmap
