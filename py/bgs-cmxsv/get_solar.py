#!/usr/bin/env python                                                                                                                                                                                                                    
from    __future__               import  absolute_import, division, print_function

import  time
import  pickle
import  numpy as np

from    pkg_resources import resource_filename
from    scipy.interpolate import interp1d

import  os
import  time
import  glob
import  argparse
import  fitsio
import  desisurvey
import  warnings
import  ephem
import  numpy                    as      np
import  astropy.units            as      u
import  pylab                    as      pl
import  matplotlib.pyplot        as      plt
import  astropy.units            as      units

from    astropy.table            import  Table, vstack
from    scipy                    import  ndimage
from    multiprocessing          import  Pool, Array
from    desisurvey.utils         import  get_location
from    astropy.time             import  Time
from    astropy.coordinates      import  SkyCoord, EarthLocation, AltAz
from    specsim.atmosphere       import  krisciunas_schaefer, Moon
from    get_sky                  import  get_sky
from    pkg_resources            import  resource_filename
from    desiutil.iers            import  freeze_iers


freeze_iers()

mayall            = get_location()

emayall           = ephem.Observer()
emayall.lon       = ephem.degrees(mayall.lon.value * np.pi / 180.)
emayall.lat       = ephem.degrees(mayall.lat.value * np.pi / 180.)
emayall.elevation = mayall.height.value

moon              = ephem.Moon()
sun               = ephem.Sun()

def airmass(zd):
    # Airmass at given zenith distance.                                                                                                                                                                                                 
    return  (1. - 0.96 * np.sin(zd * np.pi / 180.)**2.)**-0.5

def get_solar(mjd, ras, decs):
    t                 = Time(mjd, format='mjd', scale='utc')
    emayall.date      = t.iso

    ras               = np.atleast_1d(ras)
    decs              = np.atleast_1d(decs)
    
    pos               = SkyCoord(ra = ras * u.degree, dec = decs * u.degree, frame='icrs').transform_to(AltAz(obstime=t, location=mayall))
    alt               = pos.alt.degree
    az                = pos.az.degree
    zd                = 90. - alt

    sun.compute(emayall)
    moon.compute(emayall)

    sun_alt             = sun.alt * (180. / np.pi)
    sun_sep             = desisurvey.utils.separation_matrix(ras * u.deg, decs * u.deg, [sun.ra] * u.deg, [sun.dec] * u.deg)[0][0].value

    moon_alt            = moon.alt * (180. / np.pi)
    moon_frac           = moon.moon_phase
    moon_sep            = desisurvey.utils.separation_matrix(ras * u.deg, decs * u.deg, [moon.ra] * u.deg, [moon.dec] * u.deg)[0][0].value
    moon_phase          = 180. * np.arccos(2.*moon_frac - 1.) / np.pi # deg.                                                                                                                                                          \
                        
    moon_zd             = (90. - moon_alt)
    X                   = airmass(zd)

    result              = Table()
    
    result['RA']        = ras
    result['DEC']       = decs

    result['MJD']       = mjd
    
    result['AIRMASS']   = X

    result['ALT']       = alt
    result['AZ']        = az

    result['SUNALT']    = sun_alt
    result['SUNSEP']    = sun_sep

    result['MOONALT']   = moon_alt
    result['MOONSEP']   = moon_sep
    result['MOONFRAC']  = moon_frac
    result['MOONPHASE'] = moon_phase
    result['MOONZD']    = moon_zd
    
    return  result

if __name__ == '__main__':
    ras  = np.array([40., 55., 86.])
    decs = np.array([35., 25., 10.]) 
    
    for day in range(10):     
        result = get_solar(59232. + day, ras, decs)

        result.pprint()
    
