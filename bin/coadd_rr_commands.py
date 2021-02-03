#!/usr/bin/env python
from    __future__           import division, print_function

import  ephem
import  fitsio
import  healpy
import  desisurvey.config
import  sys, os, glob, time, warnings
import  numpy                as     np
import  astropy.units        as     u
import  desisurvey.utils     as     dutils

from    pathlib              import Path
from    astropy.time         import Time
from    astropy.table        import Table, vstack, hstack
from    astropy.coordinates  import Angle
from    desisurvey.utils     import get_airmass
from    desiutil             import dust
from    astropy.coordinates  import SkyCoord
from    pkg_resources        import resource_filename

# os.system('source ./env.sh')

np.random.seed(314)

'''
Generate coadds and redshifts other than nightly. 
'''

# Run date to timestamp output. 
date          = '20210202'

redux_dir     = '/global/cfs/cdirs/desi/spectro/redux/daily/'
output_dir    = '/global/cscratch1/sd/mjwilson/desi/SV1/spectra/daily/exposures/'

# number of exposures in a coadded; 1 for single-exposure coadd                                                                                                                                                                             
ALL           = False   # Overrides n_exp.

n_exp         = 1
n_node        = 1
nside         = 512

sampling      = 0.10 

overwrite     = False
archetypes    = False
verbose       = False

#
cpath         = '/global/cfs/cdirs/desi/survey/observations/SV1/sv1-exposures.fits'
cond          = Table.read(cpath)
cond          = cond[cond['TARGETS'] == 'BGS+MWS']

nights        = np.unique(cond['NIGHT'].data).astype(str)
tiles         = np.unique(cond['TILEID'].data)

npetal        = 10
petals        = list(range(npetal))

print('\n\n\t---  TILES  ---')

for tile in tiles:
  print('\t    {}'.format(tile))

print('\n\t---  NIGHTS  ---')

for night in nights:
  print('\t    {}'.format(night))

print('\n')

if overwrite | (not os.path.isfile(output_dir + '/bgs_allcframes_{}.fits'.format(date))):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
  
    ################################## Get list of exposures ##################################
    exposure_dir_list = []
    
    for obsdate in nights:
        exposure_dir_list += glob.glob(os.path.join(redux_dir, 'exposures', obsdate, '*'))
        
    ################################## Get list of all science cframes ##################################  
    cframe_list = []

    # Get a list of all science exposures.
    for exposure_dir in exposure_dir_list:      
        cframe_list_tmp = glob.glob(os.path.join(exposure_dir, 'cframe-*'))

        if len(cframe_list_tmp) > 0:
            if tiles is None:
                cframe_list += cframe_list_tmp

            else:  
                # only need to check one cframe file in the exposure
                with fitsio.FITS(cframe_list_tmp[0]) as f:
                    if f[0].read_header()['TILEID'] in tiles:
                        cframe_list += cframe_list_tmp
                    
    cframe_list           = sorted(cframe_list)

    # Gather exposure/petal information
    cframes               = Table()

    cframes['cframe']     = np.array(cframe_list)
    cframes['night']      = '                                        '
    cframes['mjd']        = np.zeros(len(cframes), dtype=np.float) - 1.0 
    cframes['lat']        = np.zeros(len(cframes), dtype=np.float)
    cframes['lon']        = np.zeros(len(cframes), dtype=np.float)
    cframes['elv']        = np.zeros(len(cframes), dtype=np.float)
    cframes['tileid']     = np.zeros(len(cframes), dtype=int)
    cframes['expid']      = np.zeros(len(cframes), dtype=int)
    cframes['exptime']    = np.zeros(len(cframes), dtype=np.float)
    cframes['camera']     = '                                        '
    cframes['program']    = '                                        '
    cframes['petal_loc']  = -1 * np.ones(len(cframes), dtype=np.int32)
    cframes['specgrph']   = np.zeros(len(cframes), dtype=np.int)
    cframes['ra']         = np.zeros(len(cframes), dtype=np.float)
    cframes['dec']        = np.zeros(len(cframes), dtype=np.float)

    for index, cframe in enumerate(cframes['cframe']):
        with fitsio.FITS(cframe) as f:
            header                         = f[0].read_header()

            if verbose:
                print(header)
           
            cframes['mjd'][index]          = header['MJD-OBS']
            cframes['night'][index]        = header['NIGHT']
            cframes['tileid'][index]       = header['TILEID']
            cframes['expid'][index]        = header['EXPID']
            cframes['camera'][index]       = header['CAMERA'].strip()[0]
            cframes['petal_loc'][index]    = int(header['CAMERA'].strip()[1])
            cframes['program'][index]      = header['PROGRAM']
            cframes['lat'][index]          = header['OBS-LAT']
            cframes['lon'][index] 	   = header['OBS-LONG']
            cframes['elv'][index] 	   = header['OBS-ELEV']
            cframes['exptime'][index]      = header['EXPTIME']
            cframes['ra'][index]           = header['SKYRA']
            cframes['dec'][index]          = header['SKYDEC']
            cframes['specgrph'][index]     = header['SPECGRPH']

    uids, cnts = np.unique(cframes['expid'], return_counts=True)

    cframes.sort(('tileid', 'petal_loc'))

    ##  cframes.pprint(max_width=-1)

    cframes.write(output_dir + '/bgs_allcframes_{}.fits'.format(date), format='fits', overwrite=True)

    exit(0)
   
## 
cframes = Table.read(output_dir + '/bgs_allcframes_{}.fits'.format(date))

##  Not in Blanc.
cframes = cframes[cframes['night'].data.astype(np.int) > 20201223]

## Remove stash.
keep = ['stash' not in x for x in cframes['cframe']]
cframes = cframes[keep]

print('Found {} nights.'.format(len(np.unique(cframes['night']))))
print('Found {} tiles.'.format(len(np.unique(cframes['tileid']))))
print('Found {} exps.'.format(len(np.unique(cframes['expid']))))
print('Found {} exptimes.'.format(np.unique(np.round(cframes['exptime'].data, decimals=0))))
print('\n\n')
print(np.unique(cframes['night'].data))

# Sanity check: each petal must have three cframe files.                                                                                                                                                                                  
for expid in np.unique(cframes['expid']):                                                                                                                                                                                                 
  mask_expid = cframes['expid'] == expid                                                                                                                                                                                                
  
  for petal_loc in petals:                                                                                                                                                                                                              
    mask = mask_expid & (cframes['petal_loc'] == petal_loc)                                                                                                                                                                           
    
    if (np.sum(mask) > 0) & (np.sum(mask) != 3):                                                                                                                                                                                      
      print('Missing arm.')
      print(cframes[mask])

## 
output_argument_list = []

if not ALL:
  if (not overwrite) and (os.path.isfile(output_dir + "/scripts/commands_spectra_nexp_{}_{}.sh".format(n_exp, night)) | os.path.isfile(output_dir + "/scripts/commands_rr_nexp_{}_{}.sh".format(n_exp, night))):
    raise  ValueError('Overwrite=True required to remove exisiting files.')

  output_file   = open(output_dir + "/scripts/commands_spectra_nexp_{}_{}.sh".format(n_exp, night), "w")
  output_rrfile = open(output_dir + "/scripts/commands_rr_nexp_{}_{}.sh".format(n_exp, night),    "w")
  
else:
  if (not overwrite) and (os.path.isfile(output_dir + "/scripts/commands_spectra_allexp_{}.sh".format(night)) | os.path.isfile(output_dir + "/scripts/commands_rr_allexp_{}.sh".format(night))):
    raise  ValueError('Overwrite=True required to remove exisiting files.')
  
  output_file   = open(output_dir + "/scripts/commands_spectra_allexp_{}.sh".format(night), "w")
  output_rrfile = open(output_dir + "/scripts/commands_rr_allexp_{}.sh".format(night),    "w")

Path(output_dir + '/scripts/').mkdir(parents=True, exist_ok=True)

##
draws = []

##
for night in nights:
  nightframes   = cframes[cframes['night'] == night]

  for tileid in np.unique(nightframes['tileid']):
    for petal_loc in petals: 
        mask                    = (cframes['tileid'] == tileid) & (cframes['petal_loc'] == petal_loc) & (cframes['night'] == night)

        # Check b camera for simplicity, only reduce if all three cameras are present.  
        mask                   &= (cframes['camera'] == 'b')

        cframe1                 = cframes[mask]

        if not ALL:
          if (np.count_nonzero(mask) < n_exp):
            print('\t\t {} exposures is not enough for TILEID {}, NIGHT {} PETAL_LOC {} for a {} exp. coadd.'.format(np.sum(mask), tileid, night, petal_loc, n_exp))
            continue

          else:
            print('\t {} exposures for TILEID {}, NIGHT {} PETAL_LOC {} (for a {} exp. coadd).'.format(np.sum(mask), tileid, night, petal_loc, n_exp))
          
          # Skip the exposures that do not make the split.
          cframe1               = cframe1[:len(cframe1) - len(cframe1) % (n_exp)]

          nsplit                = len(cframe1)//(n_exp)

          subset_split          = np.split(np.arange(len(cframe1)), nsplit)

        else:
          if (np.sum(mask) == 0):
            print('\n# No exposures for TILEID {}, PETAL_LOC {}.\n'.format(tileid, petal_loc))
            continue
          
          nsplit                = 1
          subset_split          = np.split(np.arange(len(cframe1)), nsplit)
          
        for subset_index in range(len(subset_split)):
            subset              = cframe1[subset_split[subset_index]]
            input_argument      = ''

            for index in range(len(subset)):
                exposure_dir    = os.path.dirname(subset['cframe'][index])

                # .replace(redux_dir, '$REDUXDIR/')
                input_argument += os.path.join(exposure_dir, 'cframe-[brz]{}-*.fits ').format(petal_loc)

            if (not ALL) & (n_exp == 1):
                exposure        = os.path.basename(exposure_dir)

                output_argument = os.path.join(output_dir, 'NEXP{}'.format(n_exp), str(tileid), 'spectra-{}-{}-{}.fits'.format(petal_loc, tileid, exposure))

            elif not ALL:
                output_argument = os.path.join(output_dir, 'NEXP{}'.format(n_exp), str(tileid), 'spectra-{}-{}-{}exp-subset-{}.fits'.format(petal_loc, tileid, exposure, subset_index))

            else:
                output_argument = os.path.join(output_dir, 'ALL', str(tileid), night, 'spectra-{}-{}-allexp.fits'.format(petal_loc, tileid))
                
            output_argument_list.append(output_argument)

            sample = np.random.uniform(0.0, 1.0)

            draws.append(sample)

            if sample > sampling:
              continue
            
            if os.path.isfile(output_argument) and (not overwrite):
                print('\nWarning: {} already exists!\n'.format(output_argument))
                continue

            output_file.write('desi_group_spectra --inframes {} --outfile {}\n'.format(input_argument, output_argument))

draws = np.array(draws)
            
##
output_file.close()
        
for output_argument, sample in zip(output_argument_list, draws):
    if sample > sampling:
      continue
  
    rrdesi_argument_redrock = output_argument.replace('spectra-', 'rr-').replace('.fits', '.h5')
    rrdesi_argument_zbest   = output_argument.replace('spectra-', 'zbest-')

    if os.path.isfile(rrdesi_argument_redrock) & os.path.isfile(rrdesi_argument_zbest):
        print('\nWarning: {} already exists!\n'.format(rrdesi_argument_redrock))
        continue
      
    cmd                     = 'srun -N {} -n {} -c 2 rrdesi_mpi -o {} -z {} {}'.format(n_node, 32 * n_node, rrdesi_argument_redrock, rrdesi_argument_zbest, output_argument)

    if archetypes:
      cmd                   = cmd + ' --archetypes /global/common/software/desi/cori/desiconda/20190804-1.3.0-spec/code/redrock-archetypes/master/\n'    

    else:
      pass
    
    # output_argument
    # /global/cscratch1/sd/mjwilson/desi/SV1/spectra/daily/exposures/NEXP1/80612/spectra-0-80666-00068646.fits

    # $OUTDIR/65008/redrock-2-65008-00055456.log
    
    cmd += '  &> {}'.format(rrdesi_argument_redrock.replace('.h5', '.log').replace('rr-', 'redrock-'))

    cmd += '\n'
    
    output_rrfile.write(cmd)

print('\n\n{} exps results in {} rr commands ({}% sampling) in bash script to {}/scripts.'.format(len(output_argument_list), np.count_nonzero(draws <= sampling), 100. * sampling, output_dir))    
print('\n\nDone.\n\n')
