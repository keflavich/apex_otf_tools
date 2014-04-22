import numpy as np
from pyspeckit.spectrum.readers import read_class
from astropy import coordinates
from astropy import units as u

def load_apex_cube(apex_filename='data/E-085.B-0964A-2010_merge.apex',
                   skip_data=False, DEBUG=False, downsample_factor=None):
    spectra,headers,indices = read_class.read_class(apex_filename, start=1024,
                                                    DEBUG=DEBUG,
                                                    skip_data=skip_data,
                                                    downsample_factor=downsample_factor)

    for h,i in zip(headers,indices):
        h.update(i)

    return spectra,headers,indices


def select_apex_data(spectra,headers,indices, sourcename=None,
                     shapeselect=None, tsysrange=None, rchanrange=None,
                     xtel=None,
                     skip_data=False):

    print "Determining RA/Dec"
    ra,dec = zip(*[(h['RA']+h['RAoff']/np.cos(h['DEC']/180.*np.pi),
                    h['DEC']+h['DECoff']) for h in headers])
    print "Determining Galactic coordinates"
    gal = coordinates.ICRS(np.array(ra),np.array(dec),unit=(u.deg,u.deg)).galactic
    gal.l.wrap_angle = 180*u.deg
    galOK = ((gal.l.deg > -2) &
             (gal.l.deg < 2) &
             (gal.b.deg > -2) &
             (gal.b.deg < 2))

    if not skip_data:
        print "Shaping data"
        data1 = np.array(spectra)
        shapes = np.array([d.shape for d in data1])
        if shapeselect is not None:
            OKshapes = (shapes == shapeselect).squeeze()
        elif len(np.unique(shapes)) > 1:
            raise ValueError("Inconsistent shapes.")
        else:
            OKshapes = True
    else:
        OKshapes = True
    
    if sourcename is not None:
        sourceOK = np.array([h['SOURC'].strip()==sourcename for h in headers])
    else:
        sourceOK = True

    if xtel is not None:
        xtelOK = np.array([h['XTEL'].strip()==xtel for h in headers])
    else:
        xtelOK = True


    if tsysrange is not None:
        tsys = np.array([h['TSYS'] for h in headers])
        tsysOK = (tsys>tsysrange[0]) & (tsys<tsysrange[1])
    else:
        tsysOK = True

    if rchanrange is not None:
        rchan = np.array([h['RCHAN'] if 'RCHAN' in h else np.inf for h in headers])
        rchanOK = (rchan>rchanrange[0]) & (rchan<rchanrange[1])
    else:
        rchanOK = True

    allOK = galOK & OKshapes & sourceOK & tsysOK & rchanOK & xtelOK
    if allOK.sum() == 0:
        raise ValueError("Data selection yielded empty.")

    if skip_data:
        data = None
    else:
        data = np.array(data1[allOK].tolist())

    hdrs = [h for h,K in zip(headers,allOK) if K]
    gal = gal[allOK]

    return data,hdrs,gal
