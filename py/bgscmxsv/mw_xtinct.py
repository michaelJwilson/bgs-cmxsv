def mw_xtinct(ebv, band):
    # https://www.legacysurvey.org/dr9/catalogs/
    # https://arxiv.org/pdf/1012.4804.pdf
    coeffs = {'G': 3.214, 'R': 2.165, 'i': 1.592, 'Z': 1.211, 'Y': 1.064}

    Ab     = coeffs[band] * ebv

    return  10.**(-Ab / 2.5) 


if __name__ == '__main__':
    import astropy.io.fits as fits

    
    dat = fits.open('/project/projectdirs/desi/target/catalogs/dr9m/0.44.0/targets/main/resolve/dark/targets-dark-hp-136.fits')

    # dat['TARGETS'].data['EBV'][0]                0.00983336
    # dat['TARGETS'].data['MW_TRANSMISSION_G'][0]  0.97131085
    # dat['TARGETS'].data['MW_TRANSMISSION_R'][0]  0.9805829 
    # dat['TARGETS'].data['MW_TRANSMISSION_Z'][0]  0.98909205

    for band in ['G', 'R', 'Z']:
        print(band, mw_xtinct(0.00983336, band))

        # G 0.9713108328839858
        # R 0.9805828881567109
        # Z 0.9890920710822113
        
    print('\n\nDone.\n\n')
