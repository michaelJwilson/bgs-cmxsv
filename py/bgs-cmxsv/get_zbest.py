import astropy
import numpy                         as     np
import astropy.io.fits               as     fits

from   astropy.table                 import Table, join
from   pkg_resources                 import resource_filename
from   desitarget.sv1.sv1_targetmask import desi_mask         as sv1_desi_mask
from   desitarget.sv1.sv1_targetmask import bgs_mask          as sv1_bgs_mask
from   get_solar                     import get_solar
from   desisurvey.utils              import local_noon_on_date, get_date  

def get_zbest(dpath, exps=None, lunar=False):
    # Single exposure only.
    infile   = fits.open(dpath)

    zbest    = infile['ZBEST'].data
    fmap     = infile['FIBERMAP'].data

    zbest    = Table(zbest)
    fmap     = Table(fmap)
    
    tinfo         = fmap['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'MASKBITS', 'EBV', 'PHOTSYS', 'SV1_DESI_TARGET', 'SV1_BGS_TARGET', 'DESI_TARGET', 'BGS_TARGET']
    tinfo         = astropy.table.unique(tinfo, keys='TARGETID')

    expid         = np.int(dpath.split('/')[-1].split('-')[-1].replace('.fits', ''))
    fmap['EXPID'] = expid

    keep          = ['TARGETID', 'EXPID', 'FIBERSTATUS']
    
    if lunar:        
        lunar          = get_solar(exps['MJDOBS'][exps['EXPID'] == expid][0], exps['TILERA'][exps['EXPID'] == expid][0], exps['TILEDEC'][exps['EXPID'] == expid][0])
        lunar['EXPID'] = expid

        fmap           = join(fmap, lunar, keys='EXPID', join_type='left')

        keep          += lunar.dtype.names
        keep           = list(set(keep))
        
    if exps is not None:
        # fmap[TILEID, FIBERSTATUS, EXPID]                                                                                                                                                                                                 
        fmap      = join(fmap[keep], exps['EXPID', 'EXPTIME', 'B_DEPTH', 'R_DEPTH', 'Z_DEPTH', 'GFA_TRANSPARENCY_MED', 'GFA_FWHM_ASEC_MED'], keys='EXPID', join_type='left')
        fmap      = fmap[fmap['FIBERSTATUS'] == 0]
        
    deep          = join(zbest, tinfo, keys='TARGETID', join_type='left')
    deep          = join(deep,   fmap, keys='TARGETID', join_type='left')
    
    del deep['EXPID']
    del deep['FIBERSTATUS']

    assert  len(deep) == 500
    
    # Limit to BGS.                                                                                                                                                                                                                    
    badbgs   = (deep['SV1_DESI_TARGET'] & sv1_desi_mask['BGS_ANY']) == 0
            
    # nbright  = (deep['SV1_BGS_TARGET']  & sv1_bgs_mask['BGS_BRIGHT']) == 0
    # nfaint   = (deep['SV1_BGS_TARGET']  & sv1_bgs_mask['BGS_FAINT'])  == 0

    # Limit to faint  | bright bgs only.                                                                                                                                                                                                
    # badbgs = badbgs | (nbright & nfaint)                                                                                                                                                                                             
    
    # Limit to bright bgs only.                                                                                                                                                                                                        
    # badbgs = badbgs | nbright                                                                                                                                                                                                        

    #                                                                                                                                                                                                                                  
    deep     = deep[~badbgs]
    
    # https://github.com/desihub/redrock/blob/master/py/redrock/zwarning.py                                                                                                                                                            
    deep['NODATA'] = (deep['ZWARN'] & 2**9) != 0

    # print(np.count_nonzero(badfiber), np.count_nonzero(deep['NODATA']))                                                                                                                                                              
    
    # Limit to working fibers.                                                                                                                                                                                                         
    badfiber = deep['NODATA']
    deep     = deep[~badfiber]

    return  deep


if __name__ == '__main__':
    epath   = resource_filename('bgs-cmxsv', 'dat/sv1-exposures.fits')
    exps    = Table.read(epath)

    fpath   = '/global/homes/m/mjwilson/desi/SV1/spectra/exposures/NEXP1/80612/20201218/zbest-20201218-0-00068644.fits'
    zbest   = get_zbest(fpath, exps=exps, lunar=True)
    
    zbest.pprint()

    
