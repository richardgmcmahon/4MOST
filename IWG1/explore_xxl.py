"""


Astropy guidelines
    http://docs.astropy.org/en/stable/development/codeguide.html

LSST Python style guide
    https://developer.lsst.io/coding/python_style_guide.html

PEP 8 -- Style Guide for Python Code
    https://www.python.org/dev/peps/pep-0008/


cat_xmatch is now xmatch_cat

"""
from __future__ import division, print_function

# standard library functions
import inspect
import os
import sys
import time


# 3rd party functions
import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table
from astropy.coordinates import (search_around_sky, SkyCoord,
                                 match_coordinates_sky)
import astropy.units as u
from astropy.stats import mad_std


# import private functions
sys.path.append("/home/rgm/soft/python/lib/")
from librgm.plotid import plotid
from librgm.plot_radec import plot_radec
from librgm.xmatch import xmatch_cat
from librgm.xmatch import xmatch_checkplot
from librgm.xmatch import xmatch_selfcheck
from librgm import stats

# from librgm.selfmatch_check import selfmatch_check


# read config file


#def mad_med(data, axis = None):
#    return np.median(np.absolute(data - np.median(data, axis)), axis)



def rd_dessnfields(infile=None,
                   fields = ['E1', 'E2', 'S1', 'S2',
                             'C1', 'C2', 'C3', 'X1', 'X2', 'X3'],
                   debug=False):
    """
    Read the DES SN field centres

    """
    import ConfigParser

    from astropy.coordinates import SkyCoord
    from astropy import units as u

    nfields = len(fields)

    infile = 'des_snfields.cfg'
    print('Read file:', infile)
    config.read(infile)
    print('file sections:', config.sections())

    ra = np.empty(nfields, dtype=np.float64)
    dec = np.empty(nfields, dtype=np.float64)
    print(len(ra), ra)

    ifield = -1
    for field in fields:
        ifield = ifield + 1
        radec_field = config.get('DEFAULT', 'RADEC_' + field)

        radec = SkyCoord(radec_field,  unit=(u.hour, u.deg))
        print('RA, Dec:', ifield, radec.ra.deg, radec.dec.deg)

        ra[ifield] = radec.ra.deg
        dec[ifield] = radec.dec.deg

    return ra, dec


def explore_surveys_checkplots(data, rarange=None, decrange=None,
                               colname_ra='ra', colname_dec='dec',
                               showplots=False,
                               rmax=15.0,
                               markersize=2.0,
                               plotfile_label='',
                               title='', suptitle=''):

    """


    """

    ra = data[colname_ra]
    dec = data[colname_dec]

    if rarange is None:
        rarange = [np.min(ra), np.max(ra)]
    # decrange = [np.min(dec), np.max(dec)]
    if decrange is None:
        decrange = [np.min(dec), np.max(dec)]

    if plotfile_label != '':
        plotfile_label = '_' + plotfile_label
    plotfile = 'plot_radec' + plotfile_label + '.png'
    ndata = len(ra)
    markersize = min(2.0, 1000000.0/ndata)
    plot_radec(ra, dec, plotfile=plotfile,
               markersize=markersize,
               rarange=rarange, decrange=decrange,
               suptitle=suptitle)
    plt.close()

    plotfile = 'plot_selfxneighbors' + plotfile_label + '.png'
    idx_xmatch = xmatch_selfcheck.xmatch_selfcheck(
        data=data, rmax=rmax, suptitle=suptitle,
        colnames_radec=[colname_ra, colname_dec],
        units_radec=['degree', 'degree'],
        markersize=markersize,
        plotfile=plotfile)
    plt.close()

    print('len(data):', len(data))
    print('len(idx_xmatch):', len(idx_xmatch))
    ra1 = data[colname_ra]
    dec1 = data[colname_dec]
    ra2 = data[colname_ra][idx_xmatch]
    dec2 = data[colname_dec][idx_xmatch]
    plotfile = 'plot_selfxneighbors' + plotfile_label + '_xmatch_checkplot.png'
    xmatch_checkplot.xmatch_checkplot(ra1, dec1, ra2, dec2,
                     width=rmax,
                     plotfile=plotfile,
                     suptitle=suptitle)
    plt.close()

    return idx_xmatch



# do you science here
if __name__ == "__main__":

    import time
    t0 = time.time()

    import argparse
    import ConfigParser

    # read information from config file
    config = ConfigParser.ConfigParser()
    print('__file__:', __file__)
    configfile = 'explore_cfhtls_w1.cfg'
    print('Read configfile:', configfile)
    config.read(configfile)
    print('config file sections:', config.sections())

    inpath = config.get('CFHTLS-W1', 'inpath')
    print('inpath:', inpath)

    filename_cfhtls = config.get('CFHTLS-W1', 'filename_cfhtls_w1')

    filename_des = config.get('CFHTLS-W1', 'filename_des')

    filename_vhs = config.get('CFHTLS-W1', 'filename_vhs')

    filename_gaia = config.get('CFHTLS-W1', 'filename_gaia')
    filename_galex = config.get('CFHTLS-W1', 'filename_galex')
    filename_wise = config.get('CFHTLS-W1', 'filename_wise')

    filename_sdssdr8 = config.get('CFHTLS-W1', 'filename_sdssdr8')
    filename_sdssdr9 = config.get('CFHTLS-W1', 'filename_sdssdr9')

    filename_xxl_sourcelist = config.get('CFHTLS-W1', 'filename_xxl_sourcelist')

    rarange_default = [34.0, 37.0]
    decrange_default = [-6.0, -3.0]

    table_stats = False

    colnames_radec_gaia =('RA', 'DEC')
    colnames_radec_galex =('RAJ2000', 'DEJ2000')


    # use formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # so that --help lists the defaults
    description = ''' '''
    epilog = " "
    parser =  argparse.ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    infile_xxl_sourcelist = inpath + filename_xxl_sourcelist
    xxl_sourcelist = Table.read(infile_xxl_sourcelist)
    print('Number of rows', len(xxl_sourcelist))
    xxl_sourcelist.info()
    xxl_sourcelist.info('stats')


    infile_sdssdr8 = inpath + filename_sdssdr8
    print('Reading:', infile_sdssdr8)
    sdssdr8 = Table.read(infile_sdssdr8)
    print('Number of rows', len(sdssdr8))
    sdssdr8.info()
    sdssdr8.info('stats')

    # xxl_sourcelist internal checks
    colnames_radec_xxl_sourcelist = ['RA', 'DEC']
    colname_ra = colnames_radec_xxl_sourcelist[0]
    colname_dec = colnames_radec_xxl_sourcelist[1]

    data = xxl_sourcelist
    suptitle = filename_xxl_sourcelist
    plotfile_label = 'xxl_sourcelist'
    rmax = 150.0
    ndata = len(data)
    markersize = min(2.0, 100000.0/ndata)
    idx_xmatch = explore_surveys_checkplots(
                     data, rmax=rmax,
                     colname_ra=colname_ra, colname_dec=colname_dec,
                     plotfile_label=plotfile_label,
                     suptitle=suptitle)
    # key=raw_input("Enter any key to continue: ")

    plotfile = 'plot_selfxmatch_xxl_sourcelist.png'
    idx_xmatch = xmatch_selfcheck.xmatch_selfcheck(data=data, rmax=rmax,
                     suptitle=suptitle,
                     colnames_radec=[colname_ra, colname_dec],
                     units_radec=['degree', 'degree'],
                     plotfile=plotfile)



    # key=raw_input("Enter any key to continue: ")


    # xmatch XXL to SDSSDR8
    data1 = xxl_sourcelist
    colnames_radec1 = colnames_radec_xxl_sourcelist

    colnames_radec_sdssdr8 = ['RAJ2000', 'DEJ2000']
    data2 = sdssdr8
    colnames_radec2 = colnames_radec_sdssdr8

    plotfile_label='xxl_xnn_sdssdr8'
    idx1, rsep = xmatch_cat.xmatch_cat(data1=data1, data2=data2,
                                  colnames_radec1=colnames_radec1,
                                  colnames_radec2=colnames_radec2,
                                  xmatch_rmax=None,
                                  plot_rmax=None,
                                  separations=False,
                                  plotfile_label=plotfile_label,
                                  nthneighbor=1)




    # GALEX internal checks
    infile_galex = inpath + filename_galex
    galex = Table.read(infile_galex)
    data = galex
    ndata = len(data)
    print('Number of rows', len(data))
    data.info()
    data.info('stats')

    colname_ra = 'RAJ2000'
    colname_dec = 'DEJ2000'
    suptitle = filename_galex
    plotfile_label = 'galex'
    rmax = 120.0
    markersize = min(2.0, 1000000.0/ndata)
    print('markersize:', markersize)
    showplots = True
    idx_xmatch = explore_surveys_checkplots(
                     data, rmax=rmax,
                     markersize=markersize,
                     colname_ra=colname_ra, colname_dec=colname_dec,
                     plotfile_label=plotfile_label,
                     showplots=showplots,
                     suptitle=suptitle)


    # XXL match to Gaia
    try:
       hasattr(xxl_source_list)
    except:
        print('xxl_source_list is not read in')
        infile_xxl_sourcelist = inpath + filename_xxl_sourcelist
        xxl_sourcelist = Table.read(infile_xxl_sourcelist)
        print('Number of rows', len(xxl_sourcelist))
        data = xxl_sourcelist
        xxl_sourcelist.info()
        xxl_sourcelist.info('stats')

    colname_ra = 'RA'
    colname_dec = 'DEC'
    suptitle = filename_xxl_sourcelist
    plotfile_label = 'xxl_sourcelist'
    rmax = 300.0
    markersize = min(2.0, 1000000.0/ndata)
    print('markersize:', markersize)
    showplots = True
    idx_xmatch = explore_surveys_checkplots(
                     data, rmax=rmax,
                     markersize=markersize,
                     colname_ra=colname_ra, colname_dec=colname_dec,
                     plotfile_label=plotfile_label,
                     showplots=showplots,
                     suptitle=suptitle)


    data1 = xxl_sourcelist
    colnames_radec1 =  ['RA', 'DEC']

    infile_gaia = inpath + filename_gaia
    gaia = Table.read(infile_gaia)

    data2 = gaia
    colnames_radec2 = colnames_radec_gaia

    plotfile_label='xxl_source_xnn_gaia'
    idx1, rsep = xmatch_cat.xmatch_cat(data1=data1, data2=data2,
                              colnames_radec1=colnames_radec1,
                              colnames_radec2=colnames_radec2,
                              xmatch_rmax=None,
                              plot_rmax=None,
                              separations=False,
                              plotfile_label=plotfile_label,
                              nthneighbor=1)


    sys.exit()




    # CFHTLS
    infile_cfhtls = inpath + filename_cfhtls
    cfhtls = Table.read(infile_cfhtls)
    print('Number of rows', len(cfhtls))
    cfhtls.info()

    if table_stats:
        cfhtls.info('stats')
    print('Elapsed time(secs): ',time.time() - t0)

    data = cfhtls
    ndata = len(data)
    colname_ra = 'RA'
    colname_dec = 'DEC'
    suptitle = filename_cfhtls
    plotfile_label = 'cfhls'
    decrange = decrange_default
    markersize = min(2.0, 100000.0/ndata)
    idx_xmatch = explore_surveys_checkplots(
                     data,
                     colname_ra=colname_ra, colname_dec=colname_dec,
                     decrange=decrange,
                     plotfile_label=plotfile_label,
                     suptitle=suptitle)



    # WISE
    infile_wise = inpath + filename_wise
    wise = Table.read(infile_wise)
    data = wise
    print('Number of rows', len(data))
    data.info()
    data.info('stats')
    ndata = len(data)

    colname_ra = 'RA'
    colname_dec = 'DEC'

    suptitle = filename_wise
    plotfile = 'plot_radec_wise.png'
    markersize = min(1.0, 100000.0/ndata)
    rmax = 60.0

    idx_xmatch = explore_surveys_checkplots(
                     data, rmax=rmax,
                     colname_ra=colname_ra, colname_dec=colname_dec,
                     plotfile_label=plotfile_label,
                     suptitle=suptitle)



    infile_vhs = inpath + filename_vhs
    vhs = Table.read(infile_vhs)
    data = vhs
    print('Number of rows', len(data))
    data.info()
    data.info('stats')

    colname_ra = 'RA'
    colname_dec = 'DEC'
    ra = data[colname_ra]
    dec = data[colname_dec]

    rarange= [np.min(ra), np.max(ra)]
    decrange = [np.min(dec), np.max(dec)]
    suptitle = filename_vhs
    plotfile = 'plot_radec_vhs.png'
    ndata = len(ra)
    markersize = min(1.0, 100000.0/ndata)
    plot_radec(ra, dec, plotfile=plotfile,
               markersize=markersize,
               rarange=rarange, decrange=decrange,
               suptitle=suptitle)


    plt.close()
    rmax = 60.0
    idx_xmatch = xmatch_selfcheck(
                     data=data, rmax=rmax, suptitle=suptitle,
                     colnames_radec=[colname_ra, colname_dec],
                     units_radec=['degree', 'degree'])


    ra1 = data[colname_ra]
    dec1 = data[colname_dec]
    ra2 = data[colname_ra][idx_xmatch]
    dec2 = data[colname_dec][idx_xmatch]
    plotfile = 'plot_selfxmatch_vhs_checkplot.png'
    xmatch_checkplot.xmatch_checkplot(ra1, dec1, ra2, dec2,
                     width=rmax,
                     plotfile=plotfile,
                     suptitle=suptitle)
    plt.close()


    # VHS Primary objects
    # WHERE (priOrSec = 0 OR priOrSec = frameSetID)
    infile_vhs = inpath + filename_vhs
    vhs = Table.read(infile_vhs)
    iprimary = vhs['PRIM'] == 1
    data = vhs[iprimary]
    print('Number of rows', len(data))
    data.info()
    data.info('stats')

    colname_ra = 'RA'
    colname_dec = 'DEC'
    ra = data[colname_ra]
    dec = data[colname_dec]

    rarange= [np.min(ra), np.max(ra)]
    decrange = [np.min(dec), np.max(dec)]
    suptitle = filename_vhs + ':Primary sources'
    plotfile = 'plot_radec_vhsprimary.png'
    ndata = len(ra)
    markersize = min(1.0, 100000.0/ndata)
    plot_radec(ra, dec, plotfile=plotfile,
               markersize=markersize,
               rarange=rarange, decrange=decrange,
               suptitle=suptitle)


    plt.close()
    rmax = 60.0
    selfxmatch_check(data=data, rmax=rmax, suptitle=suptitle,
                     colnames_radec=[colname_ra, colname_dec],
                     units_radec=['degree', 'degree'])



    # SDSS DR9
    infile_sdssdr9 = inpath + filename_sdssdr9
    sdssdr9 = Table.read(infile_sdssdr9)
    data = sdssdr9
    print('Number of rows', len(data))
    data.info()
    data.info('stats')

    colname_ra = 'RAJ2000'
    colname_dec = 'DEJ2000'
    ra = data[colname_ra]
    dec = data[colname_dec]

    suptitle = filename_sdssdr9
    plotfile = 'plot_radec_sdssdr9.png'
    ndata = len(ra)
    markersize = min(1.0, 100000.0/ndata)
    plot_radec(ra, dec, plotfile=plotfile,
               markersize=markersize,
               rarange=rarange, decrange=decrange,
               suptitle=suptitle)


    plt.close()
    rmax = 30.0
    selfxmatch_check(data=data, rmax=rmax, suptitle=suptitle,
                     colnames_radec=[colname_ra, colname_dec],
                     units_radec=['degree', 'degree'])


    # DES data
    infile_des = inpath + filename_des
    print('Reading:', infile_des)
    des = Table.read(infile_des)
    data = des
    des.info()
    des.info('stats')

    colname_ra = 'RA'
    colname_dec = 'DEC'
    ra = data[colname_ra]
    dec = data[colname_dec]

    rarange= [np.min(ra), np.max(ra)]
    decrange = [np.min(dec), np.max(dec)]
    suptitle = infile_des
    plotfile = 'plot_radec_des.png'
    plot_radec(ra, dec, plotfile=plotfile,
               markersize=0.1,
               rarange=rarange, decrange=decrange,
               suptitle=suptitle)


    plt.close()
    rmax = 30.0
    selfxmatch_check(data=data, rmax=rmax, suptitle=suptitle,
                     colnames_radec=[colname_ra, colname_dec],
                     units_radec=['degree', 'degree'])


    infile_xxl_sourcelist = inpath + filename_xxl_sourcelist
    xxl_sourcelist = Table.read(infile_xxl_sourcelist)
    print('Number of rows', len(xxl_sourcelist))
    data = xxl_sourcelist
    xxl_sourcelist.info()
    xxl_sourcelist.info('stats')




    infile_dr7qso = config.get('CFHTLS-W1', 'infile_dr7qso')
    dr7qso = Table.read(infile_dr7qso)
    print('Number of rows', len(dr7qso))
    dr7qso.info()
    dr7qso.info('stats')

    infile_dr12qso = config.get('CFHTLS-W1', 'infile_dr12qso')
    dr12qso = Table.read(infile_dr12qso)
    print('Number of rows', len(dr12qso))
    dr12qso.info()
    dr12qso.info('stats')


    ra = xxl_sourcelist['RA']
    dec = xxl_sourcelist['DEC']
    rarange= [np.min(ra), np.max(ra)]
    decrange = [np.min(dec), np.max(dec)]
    suptitle = filename_xxl_sourcelist
    plotfile = 'plot_radec_xxl_sourcelist_all.png'
    plot_radec(ra, dec, plotfile=plotfile,
               rarange=rarange, decrange=decrange,
               suptitle=suptitle)

    inpath_ozdes = "/home/rgm/Projects/DES/OzDES/"
    filename_ozdes = "OZDES_QSO_20160627.fits"
    infile_ozdes = inpath_ozdes + filename_ozdes
    ozdes_dr1qso = Table.read('/data/desardata/OzDES/OzDES_QSO_20160916.fits')
    ozdes_dr1qso.info()


    rarange= [33.0, 38.0]
    decrange = [-8.0, -3.0]
    itest = (ra > rarange[0]) & (ra < rarange[1]) & \
            (dec > decrange[0]) & (dec < decrange[1])
    ra = ra[itest]
    dec = dec[itest]
    plotlabel='XXL:'
    plot_radec(ra, dec, plotfile=plotfile,
               plotlabel=plotlabel,
               marker='o', color='blue', markersize=2.0,
               rarange=rarange, decrange=decrange,
               suptitle=suptitle)


    ra = ozdes_dr1qso['RA']
    dec = ozdes_dr1qso['DEC']
    plotfile = 'plot_radec_xxl_sourcelist_ozdes.png'
    itest = (ra > rarange[0]) & (ra < rarange[1]) & \
            (dec > decrange[0]) & (dec < decrange[1])
    ra = ra[itest]
    dec = dec[itest]
    title = filename_ozdes
    plotlabel='OzDES:'
    plot_radec(ra, dec,
               plotlabel=plotlabel,
               plotfile=plotfile, overplot=True,
               rarange=rarange, decrange=decrange,
               marker='o', color='red', markersize=4.0,
               title=title, suptitle=suptitle)


    ra = dr7qso['RA']
    dec = dr7qso['DEC']
    plotfile = 'plot_radec_xxl_sourcelist_ozdes_dr7qso.png'
    itest = (ra > rarange[0]) & (ra < rarange[1]) & \
            (dec > decrange[0]) & (dec < decrange[1])
    ra = ra[itest]
    dec = dec[itest]
    title = filename_ozdes
    plotlabel='dr7qso:'
    plot_radec(ra, dec, plotlabel=plotlabel,
               plotfile=plotfile, overplot=True,
               rarange=rarange, decrange=decrange,
               marker='+', color='orange', markersize=3.0,
               title=title, suptitle=suptitle)


    ra = dr12qso['RA']
    dec = dr12qso['DEC']
    plotfile = 'plot_radec_xxl_sourcelist_ozdes_dr7qso_dr12qso.png'
    itest = (ra > rarange[0]) & (ra < rarange[1]) & \
            (dec > decrange[0]) & (dec < decrange[1])
    ra = ra[itest]
    dec = dec[itest]
    title = filename_ozdes
    plotlabel='dr12qso:'
    plot_radec(ra, dec,
               plotlabel=plotlabel,
               plotfile=plotfile, overplot=True,
               rarange=rarange, decrange=decrange,
               marker='+', color='green', markersize=3.0,
               title=title, suptitle=suptitle)


    sn_fields = ['X1', 'X2', 'X3']
    ra_dessn_fields, dec_dessn_fields = rd_dessnfields(fields=sn_fields,
                                                       debug=True)
    ifield = -1
    for ra_centre in ra_dessn_fields:
        ifield = ifield + 1
        dec_centre = dec_dessn_fields[ifield]

        ra_min = ra_centre - 1.25
        ra_max = ra_centre + 1.25
        rarange = [ra_min, ra_max]

        dec_min =  dec_centre - 1.25
        dec_max =  dec_centre + 1.25
        decrange = [dec_min, dec_max]


        print('RA range:', rarange)
        print('Dec range:', decrange)
        ra = ozdes_dr1qso['RA']
        dec = ozdes_dr1qso['DEC']
        plotfile = 'plot_radec_xxl_sourcelist_ozdes_' + sn_fields[ifield] \
            + '.png'
        itest = (ra > rarange[0]) & (ra < rarange[1]) & \
                (dec > decrange[0]) & (dec < decrange[1])
        ra = ra[itest]
        dec = dec[itest]
        title = filename_ozdes
        plotlabel='OzDES_' + sn_fields[ifield] + ': '


        plot_radec(ra, dec,
                   plotlabel=plotlabel,
                   plotfile=plotfile, overplot=False,
                   rarange=rarange, decrange=decrange,
                   marker='o', color='red', markersize=4.0,
                   title=title, suptitle=suptitle)

        # Note this has no cos(dec) correction
        # Also, it does plot anything!
        # radius = 1.0
        # circle1 = plt.Circle((ra_centre, dec_centre), radius, color='red',
        #                    fill=False)
        # plt.add_patch(circle1)

        # this will allow cos(dec) correction if you loop through thetas
        # could be array based
        radius = 1.0
        theta = np.linspace(0, 2*np.pi, 100)
        x = ra_centre + (radius*np.cos(theta))
        y = dec_centre + (radius*np.sin(theta))
        plt.plot(x, y)

        radius = 1.05
        theta = np.linspace(0, 2*np.pi, 100)
        x = ra_centre + (radius*np.cos(theta))
        y = dec_centre + (radius*np.sin(theta))
        plt.plot(x, y)

        radius = 1.1
        theta = np.linspace(0, 2*np.pi, 100)
        x = ra_centre + (radius*np.cos(theta))
        y = dec_centre + (radius*np.sin(theta))
        plt.plot(x, y)

        ra = xxl_sourcelist['RA']
        dec = xxl_sourcelist['DEC']
        itest = (ra > rarange[0]) & (ra < rarange[1]) & \
                (dec > decrange[0]) & (dec < decrange[1])
        ra = ra[itest]
        dec = dec[itest]

        plot_radec(ra, dec,
                   plotlabel=plotlabel,
                   plotfile=plotfile, overplot=True,
                   rarange=rarange, decrange=decrange,
                   marker='o', color='blue', markersize=2.0,
                   title=title, suptitle=suptitle)
