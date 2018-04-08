
from __future__ import (division, print_function)

import sys
import time

import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table
from astropy.coordinates import (SkyCoord, search_around_sky,
                                 match_coordinates_sky)
from astropy import units as u

sys.path.append("/home/rgm/soft/python/lib/")
from librgm.plotid import plotid

def plot_aen(data,
             title=None, suptitle=None):
    """


    """
    itest = (data['cl'] == 3)

    xdata = data['ASTROMETRIC_EXCESS_NOISE'][itest]

    ndata = len(xdata)
    plt.hist(xdata, bins=100,
             label='SDSS:Non-stellar:' + str(ndata),
             alpha=0.3)

    # plt.show()

    itest = (data['cl'] == 6)
    xdata = data['ASTROMETRIC_EXCESS_NOISE'][itest]

    ndata = len(xdata)
    plt.hist(xdata, bins=100,
             label='SDSS:Stellar: ' + str(ndata),
             alpha=0.3)

    plt.grid()
    plt.legend()
    plt.xlabel('Excess Astrometric Noise (EAN)')

    if title is not None:
        plt.title(title)

    if suptitle is not None:
        plt.suptitle(suptitle)

    plotid()

    plt.savefig('tmp_fig1.png')

    plt.show()
    plt.close()

    itest = (data['cl'] == 3)

    xdata = data['ASTROMETRIC_EXCESS_NOISE'][itest]

    ndata = len(xdata)
    print('xdata range:', np.min(xdata), np.max(xdata))
    key=raw_input("Enter any key to continue: ")
    nbins = 100
    weights = (xdata * 0.0) + 0.2
    x, bins, patches = plt.hist(xdata,
        bins=nbins, range=(0,20.0),
        normed=True,
        label='SDSS:Non-stellar:' + str(ndata),
        alpha=0.3)


    itest = (data['cl'] == 6)
    xdata = data['ASTROMETRIC_EXCESS_NOISE'][itest]

    ndata = len(xdata)
    x, bins, p = plt.hist(xdata, bins=nbins, range=(0,20.0),
             normed=True,
             label='SDSS:Stellar:' + str(ndata),
             alpha=0.3)
    print('bin height range:', np.min(p), np.max(p))


    if title is not None:
        plt.title(title + ' (normed=True)')

    if title is None:
        plt.title('(normed=True)')

    if suptitle is not None:
        plt.suptitle(suptitle)

    plt.grid()
    plt.legend()
    plt.xlabel('Excess Astrometric Noise (EAN)')
    plotid()

    plt.savefig('plot_normalised_total.png')
    plt.show()
    plt.close()


    # normalise to peak bin
    itest = (data['cl'] == 3)

    xdata = data['ASTROMETRIC_EXCESS_NOISE'][itest]

    ndata = len(xdata)
    print('xdata range:', np.min(xdata), np.max(xdata))
    key=raw_input("Enter any key to continue: ")
    nbins = 100
    weights = (xdata * 0.0) + 0.2
    x, bins, patches = plt.hist(xdata,
        bins=nbins, range=(0,20.0),
        normed=True,
        label='SDSS:Non-stellar:' + str(ndata),
        alpha=0.3)

    print('x range:', np.min(x), np.max(x), len(x))
    print('bin range:', np.min(bins), np.max(bins), len(bins))
    print('patches range:', np.min(patches), np.max(patches), len(patches))

    # normalise the sum of heights to 1.0
    # To have the sum of height to be 1, add the following before plt.show():
    # for item in p:
    #     item.set_height(item.get_height()/sum(x))
    normhist = np.max(x)
    for item in patches:
        item.set_height(item.get_height()/normhist)

    itest = (data['cl'] == 6)
    xdata = data['ASTROMETRIC_EXCESS_NOISE'][itest]

    ndata = len(xdata)
    x, bins, patches = plt.hist(xdata, bins=nbins, range=(0,20.0),
             normed=True,
             label='SDSS:Stellar:' + str(ndata),
             alpha=0.3)

    # normalise the sum of heights to 1.0
    normhist=np.max(x)
    for item in patches:
        item.set_height(item.get_height()/normhist)

    xmax = np.max(x)
    xmax_norm = xmax / normhist

    plt.ylim(0.0, 1.1)

    if title is not None:
        plt.title(title + ' (normed bin peak = 1)')

    if title is None:
        plt.title('(normed bin peak = 1)')

    if suptitle is not None:
        plt.suptitle(suptitle)

    plt.grid()
    plt.legend()
    plt.xlabel('Excess Astrometric Noise (EAN)')
    plotid()

    plt.savefig('plot_normalised_peak.png')
    plt.show()



if __name__ == "__main__":
    """

    look at mag v mag_errs


    """
    import os, sys
    import time

    t0 = time.time()

    GALEX = True

    inpath = '/data/4most/IWG1/'
    filename = 'gaia_xmatch_sdssdr9_wise.fits'


    if GALEX:
        filename = 'gaia_xmatch_sdssdr9_wise_galex.fits'

    infile = inpath + filename

    print('Reading:', infile)
    data = Table.read(infile)
    print('Elapsed time(secs): ',time.time() - t0)

    verbose = True
    debug = False
    if verbose:
        data.info()
        data.info('stats')
        print('Elapsed time(secs): ',time.time() - t0)

    # Vizier version of sdssdr9:
    # cl = type (class) of object (3=galaxy, 6=star)

    aen = True
    if aen:
        suptitle = infile
        plot_aen(data, suptitle=suptitle)

    xdata = data['W1MPRO'] - data['W2MPRO']
    ydata = data['PHOT_G_MEAN_MAG'] - data['W1MPRO']

    ndata = len(xdata)
    plt.plot(xdata, ydata, '.',
        label=str(ndata))
    plt.suptitle(filename)
    plt.xlabel('W1MPRO - W2MPRO (WISE:Vega, WISE:Vega)')
    plt.ylabel('G - W1MPRO (Gaia:Vega, WISE:Vega)')
    plt.legend()
    plt.grid()

    plt.show()
    plt.close()


    stellar = (data['cl'] == 6)
    nonstellar = (data['cl'] == 3)

    xdata = data['W1MPRO'] - data['W2MPRO']
    ydata = data['PHOT_G_MEAN_MAG'] - data['W1MPRO']


    xdata = xdata[stellar]
    ydata = ydata[stellar]

    ndata = len(xdata)
    plt.plot(xdata, ydata, '.',
        label='Stellar:' + str(ndata))
    plt.suptitle(filename)
    plt.xlabel('W1MPRO - W2MPRO (WISE:Vega, WISE:Vega)')
    plt.ylabel('G - W1MPRO (Gaia:Vega, WISE:Vega)')

    xdata = data['W1MPRO'] - data['W2MPRO']
    ydata = data['PHOT_G_MEAN_MAG'] - data['W1MPRO']

    xdata = xdata[nonstellar]
    ydata = ydata[nonstellar]

    ndata = len(xdata)
    plt.plot(xdata, ydata, '.',
        label='Non-stellar:' + str(ndata))

    plt.legend()
    plt.grid()

    plt.show()
    plt.close()



    # GALEX
    if not GALEX:
        sys.exit()

    xdata = data['W1MPRO'] - data['W2MPRO']
    ydata = data['NUV'] - data['PHOT_G_MEAN_MAG']

    ndata = len(xdata)
    plt.plot(xdata, ydata, '.',
        label=str(ndata))
    plt.xlabel('W1MPRO - W2MPRO (WISE:Vega, WISE:Vega)')
    plt.ylabel('NUV - G  (GALEX:?. Gaia:Vega)')
    plt.legend()
    plt.grid()
    plt.show()

    plt.close()


    stellar = (data['cl'] == 6)
    nonstellar = (data['cl'] == 3)

    xdata = data['W1MPRO'] - data['W2MPRO']
    ydata = data['NUV'] - data['PHOT_G_MEAN_MAG']

    xdata = xdata[stellar]
    ydata = ydata[stellar]

    ndata = len(xdata)
    plt.plot(xdata, ydata, '.',
        label='Stellar:' + str(ndata))
    plt.suptitle(filename)
    plt.xlabel('W1MPRO - W2MPRO (WISE:Vega, WISE:Vega)')
    plt.ylabel('NUV - G  (GALEX:?. Gaia:Vega)')

    xdata = data['W1MPRO'] - data['W2MPRO']
    ydata = data['NUV'] - data['PHOT_G_MEAN_MAG']

    xdata = xdata[nonstellar]
    ydata = ydata[nonstellar]

    ndata = len(xdata)
    plt.plot(xdata, ydata, '.',
        label='Non-stellar:' + str(ndata))

    plt.legend()
    plt.grid()

    plt.show()
    plt.close()


    # Galex-WISE

    xdata = data['W1MPRO'] - data['W2MPRO']
    ydata = data['FUV'] - data['NUV']

    ndata = len(xdata)
    plt.plot(xdata, ydata, '.',
        label=str(ndata))
    plt.suptitle(filename)
    plt.xlabel('W1MPRO - W2MPRO (WISE:Vega, WISE:Vega)')
    plt.ylabel('FUV NUV  (GALEX:?. GALEX:Vega)')
    plt.legend()
    plt.grid()
    plt.show()

    plt.close()


    stellar = (data['cl'] == 6)
    nonstellar = (data['cl'] == 3)

    xdata = data['W1MPRO'] - data['W2MPRO']
    ydata = data['FUV'] - data['NUV']

    xdata = xdata[stellar]
    ydata = ydata[stellar]

    ndata = len(xdata)
    plt.plot(xdata, ydata, '.',
        label='Stellar:' + str(ndata))
    plt.suptitle(filename)
    plt.xlabel('W1MPRO - W2MPRO (WISE:Vega, WISE:Vega)')
    plt.ylabel('FUV - NUV (GALEX:?. GALEX:Vega)')

    xdata = data['W1MPRO'] - data['W2MPRO']
    ydata = data['FUV'] - data['NUV']

    xdata = xdata[nonstellar]
    ydata = ydata[nonstellar]

    ndata = len(xdata)
    plt.plot(xdata, ydata, '.',
        label='Non-stellar:' + str(ndata))

    plt.legend()
    plt.grid()

    plt.show()
    plt.close()
