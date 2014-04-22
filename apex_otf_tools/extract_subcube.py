from astropy.io import fits
from astropy import units as u
import numpy as np
from astropy import constants
from FITS_tools import cube_regrid

def extract_subcube(cubefilename, outfilename, linefreq=218.22219*u.GHz,
                    debug=False, smooth=False, vsmooth=False):
    ffile = fits.open(cubefilename)
    # I don't know why this is necessary, but there were weird bugs showing up
    # if I did not add this (some parts of the cube defaulted to 3e-319)
    ffile[0].data = ffile[0].data.astype('float32')
    hdr = ffile[0].header

    cdk = 'CD3_3' if 'CD3_3' in hdr else 'CDELT3'

    if debug:
        xarr = (np.arange(hdr['NAXIS3'])+1-hdr['CRPIX3']) * hdr[cdk] + hdr['CRVAL3']
        print xarr.min(),xarr.max()

    hdr['CTYPE3'] = 'VELO'
    cdfrq = hdr[cdk]
    crvfreq = hdr['CRVAL3']
    crpixf = hdr['CRPIX3']
    hdr[cdk] = -((cdfrq/crvfreq)*constants.c).to(u.km/u.s).value
    hdr['CRVAL3'] = 0.0
    hdr['CRPIX3'] = (linefreq.to(u.Hz).value-crvfreq)/cdfrq + crpixf
    hdr['CUNIT3'] = 'km/s'

    if debug:
        xarr = (np.arange(hdr['NAXIS3'])+1-hdr['CRPIX3']) * hdr[cdk] + hdr['CRVAL3']
        print xarr.min(),xarr.max()

    for k in ['CTYPE3','CRVAL3',cdk,'CRPIX3','CUNIT3']:
        assert ffile[0].header[k] == hdr[k]

    outhdr = hdr.copy()
    outhdr[cdk] = 1.0
    outhdr['RESTFREQ'] = linefreq.to(u.Hz).value
    outhdr['RESTFRQ'] = linefreq.to(u.Hz).value
    outhdr['CRVAL3'] = 50
    outhdr['CRPIX3'] = 199.5
    outhdr['NAXIS3'] = 400

    if smooth:
        #cubesm = gsmooth_cube(ffile[0].data, [3,2,2], use_fft=True,
        #                      psf_pad=False, fft_pad=False)
        # smoothed with 2 pixels -> sigma=10", fwhm=23"
        # this is an "optimal smooth", boosting s/n and smoothing to 36"
        # resolution.
        kw = 2 if not vsmooth else 4
        cubesm = cube_regrid.spatial_smooth_cube(ffile[0].data, kw,
                                                 use_fft=False,
                                                 numcores=4)
        cubesm = cube_regrid.spectral_smooth_cube(cubesm, 3/2.35,
                                                  use_fft=False,
                                                  numcores=4)
        ffile[0].data = cubesm

        outhdr[cdk] = 3.0
        outhdr['CRVAL3'] = 50
        outhdr['CRPIX3'] = 70
        outhdr['NAXIS3'] = 140
    
    if debug:
        xarr = (np.arange(outhdr['NAXIS3'])+1-outhdr['CRPIX3']) * outhdr[cdk] + outhdr['CRVAL3']
        print xarr.min(),xarr.max()
        xarr = (-xarr/3e5) * linefreq + linefreq
        print xarr.min(),xarr.max()
        return hdr,outhdr

    newhdu = cube_regrid.regrid_cube_hdu(ffile[0], outhdr, order=1,
                                         prefilter=False)

    newhdu.writeto(outfilename, clobber=True)
    
    return newhdu
