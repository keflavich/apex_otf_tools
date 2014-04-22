import numpy as np
from sdpy import makecube
from astropy import units as u
from astropy import constants

def add_apex_data(data,hdrs,gal, cubefilename, noisecut=np.inf, retfreq=False,
                  excludefitrange=None, varweight=False, debug=False):


    print "Data shape: ",data.shape
    if data.ndim != 2:
        raise ValueError('Data shape is NOT ok.')
    if data.shape[0] != len(hdrs):
        raise ValueError('Data and headers do not match')
    if data.shape[0] != len(gal):
        raise ValueError('Data and coords od not match')

    def data_iterator(data=data, continuum=False, fsw=False):
        shape0 = data.shape[0]
        for ii in xrange(shape0):
        #for ii in xrange(1000):
            yield data[ii,:]

    # as defined on http://www.apex-telescope.org/heterodyne/shfi/het230/lines/
    linefreq = 218222.192
    def velo_iterator(data=None, linefreq=linefreq, headers=hdrs):
        for h in headers:
            if retfreq:
                freqarr = ((np.arange(h['NCHAN'])+1-h['RCHAN']) * h['FRES'] +
                           h['FOFF'] + h['RESTF'])
                #veloarr = ((freqarr-linefreq)/linefreq * constants.c).to(u.km/u.s).value
                # needs to be in hz
                yield freqarr*1e6
            else:
                veloarr = (np.arange(h['NCHAN'])+1-h['RCHAN']) * h['VRES'] + h['VOFF']
                yield veloarr

    def coord_iterator(data=None, coordsys_out='galactic', gal=gal):
        for c in gal:
            yield c.l.deg, c.b.deg

    nhits = cubefilename+"_nhits.fits"

    makecube.add_data_to_cube(cubefilename+".fits", data=data, flatheader='header.txt',
                              cubeheader='cubeheader.txt', linefreq=218.22219, allow_smooth=True,
                              nhits=nhits,
                              data_iterator=data_iterator,
                              coord_iterator=coord_iterator,
                              velo_iterator=velo_iterator, debug=debug,
                              progressbar=True, coordsys='galactic',
                              velocity_offset=0.0, negative_mean_cut=None,
                              add_with_kernel=True, kernel_fwhm=5/3600., fsw=False,
                              diagnostic_plot_name=None, chmod=False,
                              smoothto=2,
                              noisecut=noisecut,
                              excludefitrange=None,
                              varweight=varweight,
                              continuum_prefix=None)

def make_blanks(gal, header, cubefilename, clobber=True):

    lrange = gal.l.deg.min()+15/3600.,gal.l.deg.max()+15/3600.
    brange = gal.b.deg.min()+15/3600.,gal.b.deg.max()+15/3600.
    print "Map extent: %0.2f < l < %0.2f,  %0.2f < b < %0.2f" % (lrange[0],
                                                                 lrange[1],
                                                                 brange[0],
                                                                 brange[1])

    pixsize = 5*u.arcsec
    naxis1 = (lrange[1]-lrange[0])/(pixsize.to(u.deg).value)
    naxis2 = (brange[1]-brange[0])/(pixsize.to(u.deg).value)
    restfreq = (header['RESTF']*u.MHz)
    bmaj = (1.22*10*u.m / restfreq.to(u.m,u.spectral()))*u.radian

    makecube.generate_header(np.mean(lrange), np.mean(brange), naxis1=naxis1,
                             naxis2=naxis2, naxis3=4096, coordsys='galactic',
                             ctype3='VELO-LSR',
                             bmaj=bmaj.to(u.deg).value,
                             bmin=bmaj.to(u.deg).value,
                             pixsize=pixsize.to(u.arcsec).value,
                             cunit3='km/s',
                             output_flatheader='header.txt',
                             output_cubeheader='cubeheader.txt',
                             cd3=header['VRES'],
                             crval3=-1*header['VRES']*header['RCHAN'],
                             crpix3=1,
                             clobber=True, bunit="K",
                             restfreq=restfreq.to(u.Hz).value, radio=True)

    makecube.make_blank_images(cubefilename, clobber=clobber)

def make_blanks_freq(gal, header, cubefilename, clobber=True):
    """ complete freq covg """

    lrange = gal.l.deg.min()+15/3600.,gal.l.deg.max()+15/3600.
    brange = gal.b.deg.min()+15/3600.,gal.b.deg.max()+15/3600.
    print "Map extent: %0.2f < l < %0.2f,  %0.2f < b < %0.2f" % (lrange[0],
                                                                 lrange[1],
                                                                 brange[0],
                                                                 brange[1])

    pixsize = 5*u.arcsec
    naxis1 = int((lrange[1]-lrange[0])/(pixsize.to(u.deg).value)+10)
    naxis2 = int((brange[1]-brange[0])/(pixsize.to(u.deg).value)+10)
    restfreq = (header['RESTF']*u.MHz)
    bmaj = (1.22*10*u.m / restfreq.to(u.m,u.spectral()))*u.radian
    rchan = header['RCHAN']

    #scalefactor = 1./downsample_factor
    #crpix3 = (rchan-1)*scalefactor+0.5+scalefactor/2.

    makecube.generate_header(np.mean(lrange), np.mean(brange), naxis1=naxis1,
                             naxis2=naxis2, naxis3=header['NCHAN'],
                             coordsys='galactic',
                             bmaj=bmaj.to(u.deg).value,
                             bmin=bmaj.to(u.deg).value,
                             pixsize=pixsize.to(u.arcsec).value,
                             cunit3='Hz',
                             ctype3='FREQ',
                             output_flatheader='header.txt',
                             output_cubeheader='cubeheader.txt',
                             cd3=header['FRES']*1e6,
                             crval3=restfreq.to(u.Hz).value,
                             crpix3=rchan,
                             clobber=True, bunit="K",
                             restfreq=restfreq.to(u.Hz).value, radio=True)

    makecube.make_blank_images(cubefilename, clobber=clobber)


def make_blanks_merge(cubefilename, lowhigh='low', clobber=True):
    pixsize = 5*u.arcsec
    naxis1 = 950
    naxis2 = 300
    restfreq = 218222.192*u.MHz
    bmaj = (1.22*10*u.m / restfreq.to(u.m,u.spectral()))**-1*u.radian
    cd3 = ((1*u.km/u.s)/constants.c * 218.2*u.GHz).to(u.Hz).value
    naxis3 = int(np.ceil(((1.0*u.GHz / (218.2*u.GHz) * constants.c) / (u.km/u.s)).decompose().value))

    makecube.generate_header(0.35, -0.075, naxis1=naxis1,
                             naxis2=naxis2, naxis3=naxis3, coordsys='galactic',
                             bmaj=bmaj.to(u.deg).value,
                             bmin=bmaj.to(u.deg).value,
                             pixsize=pixsize.to(u.arcsec).value,
                             cunit3='Hz',
                             ctype3='FREQ',
                             output_flatheader='header.txt',
                             output_cubeheader='cubeheader.txt',
                             cd3=cd3,
                             crval3=216.8e9 if lowhigh=='low' else 218e9,
                             crpix3=1,
                             clobber=True, bunit="K",
                             restfreq=restfreq.to(u.Hz).value, radio=True)

    makecube.make_blank_images(cubefilename, clobber=clobber)
