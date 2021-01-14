#!/usr/bin/env python
import os
import numpy as np
import astropy.io.fits as fits
import pylab as pl
import astropy
import matplotlib.pyplot as plt

from   pathlib import Path
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

write   = False

# Use truth catalogue to check z success.
truth   = True
verbose = True
qa      = True

epath   = resource_filename('bgs-cmxsv', 'dat/sv1-exposures.fits')
exps    = Table.read(epath)
exps    = exps[exps['TARGETS'] == 'BGS+MWS']

tiles   = np.unique(exps['TILEID'].data)

# TILEID  NIGHT  EXPID
conds          = exps['TILEID', 'NIGHT', 'GFA_TRANSPARENCY_MED', 'GFA_FWHM_ASEC_MED', 'GFA_SKY_MAG_AB_MED']
conds_grouped  = conds.group_by(['TILEID', 'NIGHT'])
conds_binned   = conds_grouped.groups.aggregate(np.mean)
conds          = conds_binned

exps           = exps['TILEID', 'NIGHT', 'TILERA', 'TILEDEC', 'EXPTIME', 'B_DEPTH', 'R_DEPTH', 'Z_DEPTH']
exps['NEXP']   = 1

tinfo          = astropy.table.unique(exps['TILEID', 'TILERA', 'TILEDEC'], keys='TILEID')

exps_grouped   = exps.group_by(['TILEID', 'NIGHT'])
exps_binned    = exps_grouped.groups.aggregate(np.sum)
exps           = exps_binned

exps           = exps['TILEID', 'NIGHT', 'NEXP', 'EXPTIME', 'B_DEPTH', 'R_DEPTH', 'Z_DEPTH']
exps           = join(exps, tinfo, keys='TILEID', join_type='left')

exps           = join(exps, conds, keys=['TILEID', 'NIGHT'], join_type='left')

exps.sort('TILEID')

exps['NBGSA'] = 0
exps['NBGSW'] = 0
exps['NBGSZ'] = 0

nights        = exps['NIGHT'].data.astype(np.int)
in_blanc      = nights <= 20201223

exps          = exps[in_blanc]
nights        = np.unique(exps['NIGHT'].data).astype(np.int)

# print(nights)

truth_cache = {}

def truth_test(tileid, cat, ax=None):
  if tileid not in truth_cache.keys():  
    truth       = Table.read('/global/homes/m/mjwilson/desi/SV1/spectra/truth/bgs_deep_truth_{}.fits'.format(tileid))
    truth.sort('TARGETID')

    truth_cache[tileid] = truth

  else:
    truth    = truth_cache[tileid]

  in_cat     = np.isin(truth['TARGETID'], cat['TARGETID'])
  truth      = truth[in_cat]
    
  # Both sorted by targetid.
  for x in [truth, cat]:
    assert  np.all(np.argsort(x) == np.arange(len(x)))
    
  #
  in_truth   = np.isin(cat['TARGETID'], truth['TARGETID'])

  badz       = np.zeros(len(cat)).astype(bool)

  assert  np.all(cat['TARGETID'][in_truth] == truth['TARGETID'])
  
  badz[in_truth] = np.abs(cat['Z'][in_truth] - truth['Z']) > 2. * cat['ZERR'][in_truth]
  
  if ax is not None: 
    drop = is_badz(cat[in_truth])
    
    ax.plot(truth['Z'][~drop], cat['Z'][in_truth][~drop], marker='.', c='k', lw=0.0, markersize=2)
        
    ax.set_xlabel(r'$z_{\rm{deep}}$')
    ax.set_ylabel(r'$z_{\rm{nightly}}$')
    
  return  badz

if qa:
  fig, ax = plt.subplots(1, 1, figsize=(5,5))

else:
  ax = None
  
for tileid in tiles:
  nights      = np.unique(exps[exps['TILEID'] == tileid]['NIGHT'])

  for night in nights:   
    for petal in range(10):
        # /global/cfs/cdirs/desi/users/raichoor/fiberassign-sv1/20201212/fba-080605.fits
        path         = '/global/cfs/cdirs/desi/spectro/redux/blanc/tiles/{}/{}/zbest-{}-{}-{}.fits'.format(tileid, night, petal, tileid, night)

        # print(path)
        
        if os.path.isfile(path):            
            infile   = fits.open(path)

            zbest    = infile['ZBEST'].data
            fmap     = infile['FIBERMAP'].data
            
            zbest    = Table(zbest)
            fmap     = Table(fmap)

            print('\n\n{} \t {} \t {} \t ({} \t {})'.format(tileid, night, petal, len(zbest), len(fmap)))
            
            tinfo    = fmap['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'FLUX_R', 'FIBERFLUX_R', 'PHOTSYS', 'SV1_DESI_TARGET', 'SV1_BGS_TARGET', 'DESI_TARGET', 'BGS_TARGET'] 
            tinfo    = astropy.table.unique(tinfo, keys='TARGETID')
            
            # cols = fmap.dtype.names
            # cols = [x for x in cols if 'DESI_TARGET' in x]
            # print(tileid, cols)

            deep     = join(zbest, tinfo, keys='TARGETID', join_type='left')
            
            # FIBERSTATUS cut.                                                                                                                                                                                  
            # fibstat  = fmap['TARGETID', 'FIBERSTATUS']
            # fibstat  = fibstat.group_by('TARGETID')
            # fibstat  = fibstat.groups.aggregate(np.all)

            # badfiber = fibstat['FIBERSTATUS']

            # deep     = join(deep, fibstat, keys='TARGETID', join_type='left')

            assert  len(deep) == 500

            # Limit to BGS.                                                                                                                                                                                                       
            badbgs   = (deep['SV1_DESI_TARGET'] & sv1_desi_mask['BGS_ANY']) == 0

            nbright  = (deep['SV1_BGS_TARGET'] & sv1_bgs_mask['BGS_BRIGHT']) == 0
            nfaint   = (deep['SV1_BGS_TARGET'] & sv1_bgs_mask['BGS_FAINT'])  == 0
            
            # Limit to faint | bright bgs only.
            badbgs   = badbgs | (nbright & nfaint)
            
            # Limit to bright bgs only.                                                                                                                                                                                               
            badbgs   = badbgs | nbright
            
            #
            deep     = deep[~badbgs]

            print('{} BGS ASSIGNED.'.format(len(deep)))
            
            # Assigned BGS, e.g. 3259 for 80614 here: https://data.desi.lbl.gov/desi/users/raichoor/fiberassign-sv1/sv1-per-tile/index.html#tile-nexp-design
            exps['NBGSA'][(exps['TILEID'] == tileid) & (exps['NIGHT'] == night)] += len(deep)
            
            # https://github.com/desihub/redrock/blob/master/py/redrock/zwarning.py
            # NODATA for 'malfunctioning' positioner.
            deep['NODATA'] = (deep['ZWARN'] & 2**9) != 0 
            
            # print(np.count_nonzero(badfiber), np.count_nonzero(deep['NODATA']))

            # Limit to working fibers.
            badfiber = deep['NODATA']
            deep     = deep[~badfiber]

            exps['NBGSW'][(exps['TILEID'] == tileid) & (exps['NIGHT'] == night)] += len(deep)

            print('{} BGS WORKING FIBER.'.format(len(deep)))
            
            # Good redshifts (in deep)
            badz     = is_badz(deep, verbose=verbose)
            
            if truth:
                badz_truth = truth_test(tileid, deep, ax=ax)
              
                print('LOST {} ON TRUTH CUT.'.format(np.count_nonzero(~badz & badz_truth)))
              
                badz = badz | badz_truth
  
            #
            deep['BGS_SUCCESS'] = ~badz

            if write:
                # Write only BGS (bright) on good fibers. 
                Path('/global/cscratch1/sd/mjwilson/desi/SV1/spectra/bgs-zbest/{}/{}'.format(tileid, night)).mkdir(parents=True, exist_ok=True)
            
                deep.write('/global/cscratch1/sd/mjwilson/desi/SV1/spectra/bgs-zbest/{}/{}/bgs-zbest-{}-{}-{}.fits'.format(tileid, night, petal, tileid, night), overwrite=True)
            
            deep = deep[~badz]
                        
            for label, lost in zip(['FIBER', 'BGS', 'REDSHIFT'], [badfiber, badbgs, badz]):
                print('LOST {} on {} cut'.format(np.count_nonzero(lost), label))
                
            exps['NBGSZ'][(exps['TILEID'] == tileid) & (exps['NIGHT'] == night)] += len(deep)

            print('{} GOOD BGS.'.format(len(deep)))
            
        else:
            print('Failed to retrieve {}.'.format(path)) 
            
exps['BGSSUCCESS_%'] = ['{:.2f}'.format(100. * x['NBGSZ'] / x['NBGSW']) for x in exps]

exps.sort('BGSSUCCESS_%')

print('\n\n')

exps.pprint(max_lines=-1)

if qa:
    pl.show()
 
    pl.clf()
   
    pl.semilogx(exps['B_DEPTH'], exps['BGSSUCCESS_%'].data.astype(np.float), c='b', marker='.', lw=0.0, label='B')
    pl.semilogx(exps['R_DEPTH'], exps['BGSSUCCESS_%'].data.astype(np.float), c='g', marker='.', lw=0.0, label='R')
    pl.semilogx(exps['Z_DEPTH'], exps['BGSSUCCESS_%'].data.astype(np.float), c='r', marker='.', lw=0.0, label='Z')

    pl.xlabel('EFF. DEPTH  [S]')
    pl.ylabel('BGS SUCCESS [%]')

    pl.ylim(0.0, 100.)
  
    pl.legend(frameon=False)

    pl.show()
  
print('\n\nDone.\n\n')
