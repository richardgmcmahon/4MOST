from __future__ import (division, print_function)

import sys
import time

import numpy as np

from astropy.table import Table
from astropy.coordinates import (SkyCoord, search_around_sky,
                                 match_coordinates_sky)
from astropy import units as u

import matplotlib.pyplot as plt

sys.path.append("/home/rgm/soft/python/lib/")
from librgm.plotid import plotid


def plot_radec(table=None, plotfile_prefix=None):
    """

    """
    plt.figure(num=None, figsize=(10.0, 10.0))

    xcolname = 'ra'
    ycolname = 'dec'
    xdata = table[xcolname]
    ydata = table[ycolname]

    ndata = len(table)
    xmin = np.min(xdata)
    xmax = np.max(xdata)
    ymin = np.min(ydata)
    ymax = np.max(ydata)

    print('Number of points:', ndata)
    alpha = 1.0
    markersize = 0.5
    plt.plot(xdata, ydata, '.b', ms=0.5, alpha=alpha, markeredgecolor='b',
             label=str(ndata))

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    plt.legend(loc='upper left')

    #plt.title('Gaia-DR1: COSMOS')
    plt.suptitle(infile, fontsize='medium')
    plt.xlabel(xcolname)
    plt.ylabel(ycolname)
    plotid.plotid()

    plotfile = plotfile_prefix + '_' + 'radec.png'
    print('saving :', plotfile)
    plt.savefig(plotfile)
    plt.close()


def compare_radec():
    """

    """
    # infile = '/data/desardata/Gaia/gaia_dr1_4as_dr7qso.vot'
    ra_columns = ['_RAJ2000', 'RAJ2000', 'ra_ep2000', 'ra']
    dec_columns = ['_DEJ2000', 'DEJ2000', 'dec_ep2000', 'dec']
    # infile = '/data/desardata/Gaia/dr7qso_gaiadr1_10as.vot'
    # create a dict of the RA, Dec pairs

    action = 'Comparing: _[RA, DE]J2000; [ra, dec]'
    print(action)
    ra1 = table['_RAJ2000']
    dec1 = table['_DEJ2000']
    ra2 = table['ra']
    dec2 = table['dec']
    plots = True
    results = add_columns_spherical_offsets(
        table=table, colnames=None,
        ra1=ra1, dec1=dec1, ra2=ra2, dec2=dec2,
        plots=plots, saveplots=True,
        plot_title=action,
        plot_suptitle=infile,
        plotfile_prefix=plotfile_prefix,
        verbose=True)
    print('Elapsed time(secs): ',time.time() - t0)

    sys.exit()

    print()
    action = 'Comparing: _[RA, DE]J2000, _DEJ2000; [RA, DE]J2000'
    print(action)
    ra1 = table['_RAJ2000']
    dec1 = table['_DEJ2000']
    ra2 = table['RAJ2000']
    dec2 = table['DEJ2000']
    plots = True
    results = add_columns_spherical_offsets(
        table=table, colnames=None,
        ra1=ra1, dec1=dec1, ra2=ra2, dec2=dec2,
        plots=plots, saveplots=True,
        plot_title=action,
        plot_suptitle=infile,
        plotfile_prefix=plotfile_prefix)
    print('Elapsed time(secs): ',time.time() - t0)

    print()
    action = 'Comparing: _[RA, DE]J2000; [ra, dec]_ep2000'
    print(action)
    ra1 = table['_RAJ2000']
    dec1 = table['_DEJ2000']
    ra2 = table['ra_ep2000']
    dec2 = table['dec_ep2000']
    results = add_columns_spherical_offsets(table=table,
                                            ra1=ra1, dec1=dec1,
                                            ra2=ra2, dec2=dec2,
                                            plots=plots,
                                            colnames=None,
                                            plot_title=action,
                                            plot_suptitle=infile)
    print('Elapsed time(secs): ',time.time() - t0)

    print()
    action = 'Comparing: _[RA, DE]J2000; [ra, dec]'
    print(action)
    ra1 = table['_RAJ2000']
    dec1 = table['_DEJ2000']
    ra2 = table['ra']
    dec2 = table['dec']
    results = add_columns_spherical_offsets(table=table,
                                            ra1=ra1, dec1=dec1,
                                            ra2=ra2, dec2=dec2,
                                            plots=plots,
                                            colnames=None,
                                            plot_title=action,
                                            plot_suptitle=infile)
    print('Elapsed time(secs): ',time.time() - t0)

    print()
    action = 'Comparing: [ra, dec]_ep2000; [ra, dec]'
    print(action)
    ra1 = table['ra_ep2000']
    dec1 = table['dec_ep2000']
    ra2 = table['ra']
    dec2 = table['dec']
    results = add_columns_spherical_offsets(table=table,
                                            ra1=ra1, dec1=dec1,
                                            ra2=ra2, dec2=dec2,
                                            plots=plots,
                                            colnames=None,
                                            plot_title=action,
                                            plot_suptitle=infile)
    print('Elapsed time(secs): ',time.time() - t0)

    print()


    sys.exit()



if __name__ == "__main__":
    """

    look at mag v mag_errs


    """
    import os, sys
    import time

    from astropy.table import hstack

    # could import check_matching as xm
    sys.path.append("/home/rgm/soft/python/lib/")
    from librgm.xmatch import xmatch_cat
    from librgm.xmatch import add_columns_spherical_offsets

    t0 = time.time()

    inpath = '/data/4most/IWG1/'

    filename_gaia = 'Gaia_DR1_source.fits'
    filename_galex = 'GALEX-DR5-AIS.fits'
    filename_sdssdr9 = 'sdssdr9_vizier.fits'
    filename_wise = 'AllWISE.fits'
    #filename = 'CDS_xmatch_Veron+2010_GaiaDR1_10as.vot'

    infile_gaia = inpath + filename_gaia
    infile_galex = inpath + filename_galex
    infile_sdssdr9 = inpath + filename_sdssdr9
    infile_wise = inpath + filename_wise

    print('Reading:', infile_gaia)
    gaia = Table.read(infile_gaia)
    print('Elapsed time(secs): ',time.time() - t0)

    verbose = True
    debug = False
    if verbose:
        gaia.info()
        gaia.info('stats')
        print('Elapsed time(secs): ',time.time() - t0)

    print('Reading:', infile_sdssdr9)
    sdssdr9 = Table.read(infile_sdssdr9)
    print('Elapsed time(secs): ',time.time() - t0)

    verbose = True
    debug = False
    if verbose:
        sdssdr9.info()
        sdssdr9.info('stats')
        print('Elapsed time(secs): ',time.time() - t0)

    print('Reading:', infile_wise)
    wise = Table.read(infile_wise)
    print('Elapsed time(secs): ',time.time() - t0)

    verbose = True
    debug = False
    if verbose:
        wise.info()
        wise.info('stats')
        print('Elapsed time(secs): ',time.time() - t0)


    print('Reading:', infile_galex)
    galex = Table.read(infile_galex)
    print('Elapsed time(secs): ',time.time() - t0)

    verbose = True
    debug = False
    if verbose:
        galex.info()
        galex.info('stats')
        print('Elapsed time(secs): ',time.time() - t0)


    data1 = gaia
    data2 = sdssdr9
    colnames_radec_gaia = ('RA', 'DEC')
    colnames_radec_sdssdr9 = ('RAJ2000', 'DEJ2000')
    colnames_radec1 = colnames_radec_gaia
    colnames_radec2 = colnames_radec_sdssdr9
    idxmatch, rsep = xmatch_cat(table1=data1,
                                table2=data2,
                                colnames_radec1=colnames_radec1,
                                colnames_radec2=colnames_radec2)

    ra1 =
    dec1 =
    ra2 =
    dec2 =
    add_columns_spherical_offsets(table=None,
                                  ra1=None, dec1=None,
                                  ra2=None, dec2=None,
                                  drarange=None, ddecrange=None,
                                  plots=True,
                                  colnames=None)

    print('len(gaia)', len(gaia))
    # subset of sdssdr9 sources that are matched to gaia
    gaia_xmatch_sdssdr9 = sdssdr9[idxmatch]
    print('len(gaia_xmatch_sdssdr9)', len(gaia_xmatch_sdssdr9))

    # paste Gaia columns and Gaia xmatch sdssdr9 columns together
    result = hstack([gaia, gaia_xmatch_sdssdr9])

    # create table xmatch sep < 1.0"
    itest = rsep < 1.0
    result = result[itest]
    print('len(result)', len(result))

    gaia_xmatch_sdssdr9 = result

    # Now add the WISE data
    data1 = gaia_xmatch_sdssdr9
    data2 = wise
    colnames_radec_gaia = ('RA', 'DEC')
    colnames_radec_wise = ('RA', 'DEC')
    colnames_radec1 = colnames_radec_gaia
    colnames_radec2 = colnames_radec_wise
    idxmatch, rsep = xmatch_cat.xmatch_cat(data1=data1, data2=data2,
                                           colnames_radec1=colnames_radec1,
                                           colnames_radec2=colnames_radec2)

    print('len(gaia_xmatch_sdssdr9)', len(gaia_xmatch_sdssdr9))
    # subset of sdssdr9 sources that are matched to gaia
    gaia_xmatch_wise = wise[idxmatch]
    print('len(gaia_xmatch_wise)', len(gaia_xmatch_wise))

    # paste Veron columns and ALMA source columns into array
    result = hstack([gaia_xmatch_sdssdr9, gaia_xmatch_wise])

    # create table xmatch sep < 2.0"
    itest = rsep < 2.0
    result = result[itest]
    print('len(result)', len(result))

    result.write('gaia_xmatch_sdssdr9_wise.fits', overwrite=True)

    gaia_xmatch_sdssdr9_wise = result

    data1 = gaia_xmatch_sdssdr9_wise
    data2 = galex
    colnames_radec_gaia = ('RA_1', 'DEC_1')
    colnames_radec_galex = ('RAJ2000', 'DEJ2000')
    colnames_radec1 = colnames_radec_gaia
    colnames_radec2 = colnames_radec_galex
    idxmatch, rsep = xmatch_cat.xmatch_cat(data1=data1, data2=data2,
                                           colnames_radec1=colnames_radec1,
                                           colnames_radec2=colnames_radec2)

    print('len(gaia_xmatch_sdssdr9_wise)', len(gaia_xmatch_sdssdr9_wise))
    # subset of sdssdr9 sources that are matched to gaia
    gaia_xmatch_galex= galex[idxmatch]
    print('len(gaia_xmatch_galex)', len(gaia_xmatch_galex))

    # paste Veron columns and ALMA source columns into array
    result = hstack([gaia_xmatch_sdssdr9_wise, gaia_xmatch_galex])

    # create table xmatch sep < rlimit
    itest = rsep < 3.0
    result = result[itest]
    print('len(result)', len(result))

    result.write('gaia_xmatch_sdssdr9_wise_galex.fits', overwrite=True)

    sys.exit()

    plot_radec(table=table, plotfile_prefix=plotfile_prefix)


    plt.figure(num=None, figsize=(10.0, 10.0))

    xcolname = 'phot_g_mean_mag'
    ycolname = 'phot_g_mean_flux'
    xdata = table[xcolname]
    ydata = np.log10(table[ycolname] / table[ycolname + '_error'])
    ycolname = 'log10(phot_g_mean_flux/phot_g_mean_flux_error)'

    ndata = len(table)
    xmin = np.min(xdata)
    xmax = np.max(xdata)
    ymin = np.min(ydata)
    ymax = np.max(ydata)

    print('Number of points:', ndata)
    alpha = 1.0
    markersize = 0.5
    plt.plot(xdata, ydata, '.b', ms=0.5, alpha=alpha, markeredgecolor='b',
             label=str(ndata))

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    plt.legend(loc='upper left')

    #plt.title('Gaia-DR1: COSMOS')
    plt.suptitle(infile, fontsize='small')
    plt.xlabel(xcolname)
    plt.ylabel(ycolname)
    plotid.plotid()

    plotfile = plotfile_prefix + '_' + 'mag_magerr.png'
    plt.savefig(plotfile)
    plt.close()


    # print(os.path.basename(infile))
    # print(os.path.splitext(__file__)[0])
    nth_nn = 2
    plotfile = plotfile_prefix + '_' + str(nth_nn) + 'nn_check_matches.png'
    result = xm.check_matches(infile, ['ra','dec','phot_g_mean_mag'],
                              nth_nn,
                              upperlimits=[2.0,5.0],
                              save=plotfile)

    sys.exit()

    plt.figure(num=None, figsize=(10.0, 10.0))

    xdata = table['z']
    ydata = table['rmag'] - table['phot_g_mean_mag']

    ndata = len(table)
    print('Number of points:', ndata)
    plt.plot(xdata, ydata, '.b', markeredgecolor='b', label=str(ndata))

    plt.xlim(0.0, 6.0)
    plt.ylim(-3.0, 3.0)

    plt.legend(loc='upper left')

    plt.title('Gaia-DR1 -v- DR7QSO')
    plt.xlabel('Redshift')
    plt.ylabel('r[AB] - Gaia[Vega]')
    plotid.plotid()

    plt.savefig('GaiaDR1-redshift-rG' + '.png')
    plt.close()

    plt.figure(num=None, figsize=(10.0, 10.0))

    xdata = table['z']
    ydata = (0.5*(table['gmag'] + table['rmag'])) - table['phot_g_mean_mag']

    ndata = len(table)
    print('Number of points:', ndata)
    plt.plot(xdata, ydata, '.b', markeredgecolor='b', label=str(ndata))

    plt.xlim(0.0, 6.0)
    plt.ylim(-3.0, 3.0)

    plt.legend(loc='upper left')

    plt.title('Gaia-DR1 -v- DR7QSO')
    plt.xlabel('Redshift')
    plt.ylabel('<gr>[AB] - Gaia[Vega]')
    plotid.plotid()

    plt.savefig('GaiaDR1-redshift-grG' + '.png')


    plt.figure(num=None, figsize=(10.0, 10.0))

    xdata = table['errHalfMaj']
    ydata = table['errHalfMin']

    ndata = len(table)
    print('Number of points:', ndata)
    plt.plot(xdata, ydata, '.b', markeredgecolor='b', label=str(ndata))

    # plt.xlim(0.0, 6.0)
    # plt.ylim(-3.0, 3.0)

    plt.legend(loc='upper left')

    plt.title('Gaia RA, Dec errors')
    plt.xlabel('errHalfMaj')
    plt.ylabel('errHalfMin')
    plotid.plotid()

    plt.savefig('GaiaDR1-errEllipse1' + '.png')

    plt.figure(num=None, figsize=(10.0, 10.0))

    xdata = table['errHalfMaj']/table['errHalfMin']
    ydata = table['errPosAng']

    ndata = len(table)
    print('Number of points:', ndata)
    plt.plot(xdata, ydata, '.b', markeredgecolor='b', label=str(ndata))

    # plt.xlim(0.0, 6.0)
    # plt.ylim(-3.0, 3.0)

    plt.legend(loc='upper left')

    plt.title('Gaia RA, Dec errors')
    plt.xlabel('errHalfMaj/Min')
    plt.ylabel('errPosAng')
    plotid.plotid()

    plt.savefig('GaiaDR1-errEllipse2' + '.png')
