import numpy as np


def is_badz(cat, dX2_lim=40., verbose=False):
    badz     = np.zeros(len(cat)).astype(bool) 

    labels   = ['ZWARN', 'DELTACHI2', 'SPECTYPE', 'ZRANGE', 'ZERR']
    cuts     = [cat['ZWARN'] > 0, cat['DELTACHI2'] < dX2_lim, cat['SPECTYPE'] == 'STAR', (cat['Z'] < 0.0) | (cat['Z'] > 0.6), cat['ZERR'] > (0.0005 * (1. + cat['Z']))]

    for label, cut in zip(labels, cuts):
        badz = badz | cut

        if verbose:
            print('LOST {} TO {} CUT.'.format(np.count_nonzero(cut), label))

    return  badz
