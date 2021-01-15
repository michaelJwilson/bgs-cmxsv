#!/usr/bin/env python
import os
import numpy as np
import astropy.io.fits as fits
import pylab as pl
import astropy

from   astropy.table import Table, vstack, join 

# https://github.com/desihub/desitarget/blob/master/py/desitarget/sv1/data/sv1_targetmask.yaml
from   desitarget.sv1.sv1_targetmask import desi_mask as sv1_desi_mask
from   desitarget.sv1.sv1_targetmask import bgs_mask as sv1_bgs_mask
from   pkg_resources                 import resource_filename
from   badz                          import is_badz

# [BGS_FAINT,           0, "BGS faint targets",              {obsconditions: BRIGHT|GRAY|DARK}]
# [BGS_BRIGHT,          1, "BGS bright targets",             {obsconditions: BRIGHT}]
# [BGS_FAINT_EXT,       2, "BGS faint extended targets",     {obsconditions: BRIGHT}]
# [BGS_LOWQ,            3, "BGS low quality targets",        {obsconditions: BRIGHT}]
# [BGS_FIBMAG,          4, "BGS fiber magnitude targets",    {obsconditions: BRIGHT}]

# https://github.com/desihub/desitarget/blob/master/doc/nb/target-selection-bits-and-bitmasks.ipynb

write = False

# Get SV1 exposures. 
epath = resource_filename('bgs-cmxsv', 'dat/sv1-exposures.fits')
aexps = Table.read(epath)

# Limit to BGS. 
exps  = aexps[aexps['TARGETS'] == 'BGS+MWS']

tiles = np.unique(exps['TILEID'].data)

# Limit to exposures that made it into deep.  
dexps = resource_filename('bgs-cmxsv', 'dat/blanc_deep_explist.dat') 
dexps = np.loadtxt(dexps)
dexps = exps[np.isin(exps['EXPID'], dexps[:,2])]

# Calculate summed depth for deep. 
dexps          = dexps['TILEID', 'TILERA', 'TILEDEC', 'EXPTIME', 'B_DEPTH', 'R_DEPTH', 'Z_DEPTH']

tinfo          = astropy.table.unique(dexps['TILEID', 'TILERA', 'TILEDEC'], keys='TILEID')

dexps_grouped  = dexps.group_by('TILEID')
dexps_binned   = dexps_grouped.groups.aggregate(np.sum)
dexps          = dexps_binned

dexps          = dexps['TILEID', 'EXPTIME', 'B_DEPTH', 'R_DEPTH', 'Z_DEPTH']
dexps          = join(dexps, tinfo, keys='TILEID', join_type='left')

dexps.sort('EXPTIME')

# Collect targets assigned, on a working fiber, with a good z. 
dexps['NBGSA'] = 0
dexps['NBGSW'] = 0
dexps['NBGSZ'] = 0

deep_cache     = {}

for tileid in tiles:
    deep_cache[tileid] = Table()
    
    for petal in range(10):
        # /global/cfs/cdirs/desi/users/raichoor/fiberassign-sv1/20201212/fba-080605.fits
        dpath        = '/global/cfs/cdirs/desi/spectro/redux/blanc/tiles/{}/deep/zbest-{}-{}-deep.fits'.format(tileid, petal, tileid)
    
        if os.path.isfile(dpath):            
            infile   = fits.open(dpath)

            zbest    = infile['ZBEST'].data
            fmap     = infile['FIBERMAP'].data
            
            zbest    = Table(zbest)
            fmap     = Table(fmap)
            
            print('\n\n{} \t {} \t ({} \t {})'.format(tileid, petal, len(zbest), len(fmap)))
            
            tinfo    = fmap['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'PHOTSYS', 'SV1_DESI_TARGET', 'SV1_BGS_TARGET', 'DESI_TARGET', 'BGS_TARGET'] 
            tinfo    = astropy.table.unique(tinfo, keys='TARGETID')
            
            # fmap[TILEID, FIBERSTATUS, EXPID]
            fmap     = join(fmap['TARGETID', 'EXPID', 'FIBERSTATUS'], aexps['EXPID', 'B_DEPTH', 'R_DEPTH', 'Z_DEPTH'], keys='EXPID', join_type='left')
            fmap     = fmap[fmap['FIBERSTATUS'] == 0]

            fmap_grouped = fmap.group_by('TARGETID')
            fmap         = fmap_grouped.groups.aggregate(np.sum)            
 
            deep         = join(zbest, tinfo, keys='TARGETID', join_type='left')
            deep         = join(deep,   fmap, keys='TARGETID', join_type='left')

            del deep['EXPID']
            del deep['FIBERSTATUS']
            
            assert  len(deep) == 500

            # Limit to BGS.                                                                                                                                                                                                                  
            badbgs   = (deep['SV1_DESI_TARGET'] & sv1_desi_mask['BGS_ANY']) == 0

            nbright  = (deep['SV1_BGS_TARGET']  & sv1_bgs_mask['BGS_BRIGHT']) == 0
            nfaint   = (deep['SV1_BGS_TARGET']  & sv1_bgs_mask['BGS_FAINT'])  == 0
            
            # Limit to faint | bright bgs only.
            # badbgs = badbgs | (nbright & nfaint)
            
            # Limit to bright bgs only.                                                                                                                                                                                               
            # badbgs = badbgs | nbright
            
            #
            deep     = deep[~badbgs]
            
            # Assigned BGS, e.g. 3259 for 80614 here: https://data.desi.lbl.gov/desi/users/raichoor/fiberassign-sv1/sv1-per-tile/index.html#tile-nexp-design
            dexps['NBGSA'][dexps['TILEID'] == tileid] += len(deep)
            
            # https://github.com/desihub/redrock/blob/master/py/redrock/zwarning.py
            deep['NODATA'] = (deep['ZWARN'] & 2**9) != 0 
            
            # print(np.count_nonzero(badfiber), np.count_nonzero(deep['NODATA']))

            # Limit to working fibers.
            badfiber = deep['NODATA']
            deep     = deep[~badfiber]

            del deep['NODATA']
            
            # Assigned to a working fiber.
            dexps['NBGSW'][dexps['TILEID'] == tileid] += len(deep)
                        
            badz     = is_badz(deep, dX2_lim=100.)

            # Good redshifts (in deep)  
            deep     = deep[~badz]
            
            for label, lost in zip(['FIBER', 'BGS', 'Z'], [badfiber, badbgs, badz]):
                print('LOST {} on {} cut'.format(np.count_nonzero(lost), label))
                
            # Limit to only BGS (bright) on a working fiber, with a good redshift. 
            deep_cache[tileid] = vstack((deep_cache[tileid], deep))

            dexps['NBGSZ'][dexps['TILEID'] == tileid] += len(deep)
            
    if len(deep_cache[tileid]) == 0:
        del deep_cache[tileid]

if write:
    print('\n\n')
    
    for tileid in deep_cache.keys():
        print('Writing BGS {} deep zbest.'.format(tileid))
        
        deep_cache[tileid].sort('TARGETID')
        deep_cache[tileid].write('/global/cscratch1/sd/mjwilson/desi/SV1/spectra/truth/bgs_deep_truth_{:d}.fits'.format(tileid), overwrite=True)

#
dexps['BGSSUCCESS_%'] = ['{:.2f}'.format(100. * x['NBGSZ'] / x['NBGSW']) for x in dexps]

dexps.write('/global/cscratch1/sd/mjwilson/desi/SV1/spectra/truth/summary.fits', overwrite=True)

print('\n\n')
print(dexps)
print('\n\nDone.\n\n')

