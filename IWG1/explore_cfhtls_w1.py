"""


Astropy guidelines
    http://docs.astropy.org/en/stable/development/codeguide.html

LSST Python style guide
    https://developer.lsst.io/coding/python_style_guide.html

PEP 8 -- Style Guide for Python Code
    https://www.python.org/dev/peps/pep-0008/


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
from librgm.lineno import lineno
from librgm.plot_radec import plot_radec
# from librgm import stats
from librgm.des.rd_dessnfields import rd_dessnfields

from librgm.xmatch import xmatch_cat
from librgm.xmatch import xmatch_checkplot
from librgm.xmatch import xmatch_selfcheck
from librgm.xmatch import xmatch_checkplot0

# should move to astropy mad_std
def mad_med(data, axis = None):
    return np.median(np.absolute(data - np.median(data, axis)), axis)



def explore_surveys_checkplots(data, rarange=None, decrange=None,
                               colname_ra='ra', colname_dec='dec',
                               showplots=False,
                               rmax=15.0,
                               color='black',
                               markersize=2.0,
                               plotfile_suffix='',
                               title='',
                               suptitle=''):

    """


    """

    ra = data[colname_ra]
    dec = data[colname_dec]

    if rarange is None:
        rarange = [np.min(ra), np.max(ra)]
    # decrange = [np.min(dec), np.max(dec)]
    if decrange is None:
        decrange = [np.min(dec), np.max(dec)]

    if plotfile_suffix != '':
        plotfile_suffix = '_' + plotfile_suffix
    plotfile = 'plot_radec' + plotfile_suffix + '.png'
    ndata = len(ra)
    markersize = min(2.0, 1000000.0/ndata)
    plot_radec(ra, dec, plotfile=plotfile,
               markersize=markersize,
               color=color,
               rarange=rarange, decrange=decrange,
               suptitle=suptitle)
    plt.close()

    plotfile = 'plot_selfxneighbors' + plotfile_suffix + '.png'
    idx_xmatch = xmatch_selfcheck(
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
    plotfile = 'plot_selfxneighbors' + plotfile_suffix + '_xmatch_checkplot.png'
    xmatch_checkplot(ra1, dec1, ra2, dec2,
                     width=rmax,
                     plotfile=plotfile,
                     suptitle=suptitle)
    plt.close()

    return idx_xmatch

def explore_data(data,
                 colnames_radec=None,
                 units_radec=['degree', 'degree'],
                 rarange=None,
                 decrange=None,
                 color=None,
                 infile=None,
                 filename=None,
                 title=None,
                 suptitle=None,
                 plotfile_suffix=None):

    data.info()
    print('colnames_radec:', colnames_radec)
    if 'filename' in data.meta:
        print('filename:', data.meta['filename'])

    colname_ra = colnames_radec[0]
    colname_dec = colnames_radec[1]

    ra = data[colname_ra]
    dec = data[colname_dec]

    if rarange is None:
        rarange= [np.min(ra), np.max(ra)]
    if decrange is None:
        decrange = [np.min(dec), np.max(dec)]

    if title is None:
        title = filename

    if plotfile_suffix is None:
        plotfile = 'plot_radec.png'

    if plotfile_suffix is not None:
        plotfile = 'plot_radec_' + plotfile_suffix + '.png'

    ndata = len(ra)
    markersize = min(2.0, 100000.0/ndata)
    plot_radec(ra, dec, plotfile=plotfile,
               markersize=markersize,
               rarange=rarange, decrange=decrange,
               title=title,
               color=color,
               suptitle=suptitle)


    plt.close()

    return


def plot_ozdes_xmm3dr6(xmm=None, xxl=None, ozdes_dr1=None,
    sn_fields = ['X1', 'X2', 'X3'],
    plotfile_prefix=None, suptitle='', showplots = False,
    filename_ozdes=None):
    """

    """

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

        print('RA centre, range:', ra_centre, rarange)
        print('Dec centre, range:', dec_centre, decrange)

        if ozdes_dr1 is not None:
            ra = ozdes_dr1qso['RA']
            dec = ozdes_dr1qso['DEC']

        if xmm is not None:
            ra = xmm['RAJ2000']
            dec = xmm['DEJ2000']

        plotfile = 'plot_radec_' + plotfile_suffix + sn_fields[ifield] \
                + '.png'
        itest = (ra > rarange[0]) & (ra < rarange[1]) & \
                    (dec > decrange[0]) & (dec < decrange[1])
        ra = ra[itest]
        dec = dec[itest]

        ndata = len(ra)
        print('Number of sources in RA, Dec range:', ndata)
        title = '3xmmdr6'
        if filename_ozdes is not None:
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

        plt.savefig(plotfile)
        plt.close()


        # plt.show()
        # key=raw_input("Enter any key to continue: ")

        if xxl is not None:
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
            plt.close()

    return


# do you science here
if __name__ == "__main__":

    import time
    t0 = time.time()

    import argparse
    import ConfigParser

    # read information from config file
    config = ConfigParser.ConfigParser()
    print('__file__:', __file__)
    configfile = os.path.splitext(__file__)[0] + '.cfg'
    print('Read configfile:', configfile)
    config.read(configfile)
    print('config file sections:', config.sections())

    inpath = config.get('CFHTLS-W1', 'inpath')
    print('inpath:', inpath)

    filename_cfhtls = config.get('CFHTLS-W1', 'filename_cfhtls_w1')

    filename_des = config.get('CFHTLS-W1', 'filename_des')

    filename_gaia = config.get('CFHTLS-W1', 'filename_gaia')
    filename_galex = config.get('CFHTLS-W1', 'filename_galex')
    filename_wise = config.get('CFHTLS-W1', 'filename_wise')

    filename_sdssdr8 = config.get('CFHTLS-W1', 'filename_sdssdr8')
    filename_sdssdr9 = config.get('CFHTLS-W1', 'filename_sdssdr9')

    filename_video = config.get('CFHTLS-W1', 'filename_video')
    filename_viking = config.get('CFHTLS-W1', 'filename_viking')
    filename_vhs = config.get('CFHTLS-W1', 'filename_vhs')

    filename_xxl_redshifts = config.get('CFHTLS-W1', 'filename_xxl_redshifts')
    filename_xxl_sourcelist = config.get('CFHTLS-W1', 'filename_xxl_sourcelist')

    filename_xmm3dr6 = config.get('CFHTLS-W1', 'filename_xmm3dr6')

    rarange_default = [34.0, 37.0]
    decrange_default = [-6.0, -3.0]

    colnames_radec_cfhlts = ('RA', 'DEC')
    waveband_names_cfhtls =['U', 'G', 'R', 'I', 'Z']
    magtype_cfhtls = 'MAG'

    colnames_radec_des = ('RA', 'DEC')

    colnames_radec_gaia = ('RA', 'DEC')

    colnames_radec_galex =('RAJ2000', 'DEJ2000')

    colnames_radec_vhs = ('RA', 'DEC')
    units_radec_vhs = 'deg'
    waveband_names_vhs =['J', 'H', 'K']

    colnames_radec_viking = ('RA', 'DEC')
    units_radec_viking = 'rad'

    colnames_radec_video = ('RA', 'DEC')
    units_radec_video = 'rad'

    colnames_radec_xmm3dr6 = ('RAJ2000', 'DEJ2000')

    table_stats = False

    # use formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # so that --help lists the defaults
    description = ''' '''
    epilog = " "
    parser =  argparse.ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)



    # DES
    infile_des = inpath + filename_des
    des = Table.read(infile_des)
    des.meta['filename'] = infile_des
    data = des
    print('Number of rows', len(data))
    data.info()
    data.info('stats')

    colnames_radec = colnames_radec_des

    print('RA limits:', np.min(data[colnames_radec[0]]),
          np.max(data[colnames_radec[0]]))
    print('Dec limits:', np.min(data[colnames_radec[1]]),
          np.max(data[colnames_radec[1]]))
    itest = (data[colnames_radec[0]] > rarange_default[0]) &  \
             (data[colnames_radec[0]] < rarange_default[1]) &  \
             (data[colnames_radec[1]] > decrange_default[0]) & \
             (data[colnames_radec[1]] < decrange_default[1])

    ndata = len(data)
    ndata_filtered = len(data[itest])
    print('len(data):', ndata)
    print('len(data[itest]):', ndata_filtered)

    data = data[itest]

    filename = filename_des
    plotfile_suffix=filename_des

    title = filename
    suptitle = str(ndata) + ':' + str(ndata_filtered)
    print(lineno())
    explore_data(data,
                 rarange=rarange_default,
                 decrange=decrange_default,
                 suptitle=suptitle,
                 title=title,
                 color='blue',
                 filename=filename,
                 colnames_radec=colnames_radec,
                 plotfile_suffix=plotfile_suffix)


    # 3XMM DR6 analysis
    infile_xmm3dr6 = inpath + filename_xmm3dr6

    xmm = Table.read(infile_xmm3dr6)
    xmm.meta['filename'] = infile_xmm3dr6
    data = xmm
    print('Number of rows', len(data))
    data.info()
    data.info('stats')

    plotfile_suffix='xmm3dr6_'

    plot_ozdes_xmm3dr6(xmm=xmm, sn_fields = ['E1', 'E2'],
        plotfile_prefix='xmm3dr6')
    plt.close()

    plot_ozdes_xmm3dr6(xmm=xmm, sn_fields = ['S1', 'S2'],
        plotfile_prefix='xmm3dr6')
    plt.close()

    plot_ozdes_xmm3dr6(xmm=xmm, sn_fields = ['X1', 'X2', 'X3'],
        plotfile_prefix='xmm3dr6')
    plt.close()

    plot_ozdes_xmm3dr6(xmm=xmm, sn_fields = ['C1', 'C2', 'C3'],
        plotfile_prefix='xmm3dr6')
    plt.close()

    # sys.exit()


    # CFHTLS
    infile_cfhtls = inpath + filename_cfhtls
    cfhtls = Table.read(infile_cfhtls)
    cfhtls.meta['filename'] = infile_cfhtls
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
    plotfile_suffix = 'cfhls'
    decrange = decrange_default
    markersize = min(2.0, 100000.0/ndata)
    idx_xmatch = explore_surveys_checkplots(
                     data,
                     colname_ra=colname_ra, colname_dec=colname_dec,
                     decrange=decrange,
                     color='blue',
                     plotfile_suffix=plotfile_suffix,
                     suptitle=suptitle)

    plt.close()
    colnames_radec = colnames_radec_cfhlts
    waveband_names = waveband_names_cfhtls
    magtype = magtype_cfhtls
    filename = filename_cfhtls

    print(lineno())
    key=raw_input("Enter any key to continue: ")


    # VHS: the data comes from Sergey Koposov's database and has some different
    # column names from VSA and ESO
    infile_vhs = inpath + filename_vhs
    vhs = Table.read(infile_vhs)
    vhs.meta['filename'] = infile_vhs
    data = vhs
    print('Number of rows', len(data))
    data.info()
    data.info('stats')

    plotfile_suffix='vhs'
    explore_data(data,
                 rarange=rarange_default,
                 decrange=decrange_default,
                 suptitle=None,
                 color='red',
                 filename=filename_vhs,
                 colnames_radec=colnames_radec_vhs,
                 plotfile_suffix=plotfile_suffix)

    # VHS Primary objects
    # WHERE (priOrSec = 0 OR priOrSec = frameSetID)
    infile_vhs = inpath + filename_vhs
    vhs = Table.read(infile_vhs)
    iprimary = vhs['PRIM'] == 1
    vhs = vhs[iprimary]

    data = vhs
    print('Number of rows', len(data))
    data.info()
    data.info('stats')

    suptitle = filename_vhs + ':Primary sources'
    plotfile_suffix = 'vhsprimary'
    explore_data(data,
                 rarange=rarange_default,
                 decrange=decrange_default,
                 suptitle=suptitle,
                 color='red',
                 filename=filename_vhs,
                 colnames_radec=colnames_radec_vhs,
                 plotfile_suffix=plotfile_suffix)


    # magnitude columns: ZAPERCORMAG3
    colnames_radec = colnames_radec_vhs
    waveband_names_vhs =['J', 'H', 'K']
    waveband_names = waveband_names_vhs
    magtype_vhs = 'APERCORMAG3'
    magtype = magtype_vhs
    filename = filename_vhs
    # cycle through wavebands
    for waveband in waveband_names:
        itest = (vhs[waveband + magtype] > 5.0)
        data = vhs[itest]
        plotfile_suffix = 'vhs' + '_' + waveband
        suptitle = plotfile_suffix
        explore_data(data,
                     rarange=rarange_default,
                     decrange=decrange_default,
                     suptitle=suptitle,
                     color='red',
                     filename=filename,
                     colnames_radec=colnames_radec,
                     plotfile_suffix=plotfile_suffix)

    # key=raw_input("Enter any key to continue: ")

    # VIDEO
    infile_video = inpath + filename_video
    video = Table.read(infile_video)
    video.meta['filename'] = infile_video
    colnames_radec = colnames_radec_video
    if units_radec_video == 'rad':
        video[colnames_radec[0]] = np.rad2deg(video[colnames_radec[0]])
        video[colnames_radec[1]] = np.rad2deg(video[colnames_radec[1]])
    data = video
    print('Number of rows', len(data))
    data.info()
    data.info('stats')

    filename = filename_video
    plotfile_suffix='video'
    explore_data(data,
                 rarange=rarange_default,
                 decrange=decrange_default,
                 suptitle=None,
                 color='red',
                 filename=filename,
                 colnames_radec=colnames_radec,
                 plotfile_suffix=plotfile_suffix)

    # magnitude columns: ZAPERMAG3
    waveband_names_video =['Z', 'Y', 'J', 'H', 'KS']
    magtype_video = 'APERMAG3'
    magtype = magtype_video
    # cycle through wavebands
    for waveband in waveband_names_video:
        itest = (video[waveband + magtype] > 0)
        data = video[itest]
        plotfile_suffix = 'video' + '_' + waveband
        suptitle = plotfile_suffix
        explore_data(data,
                     rarange=rarange_default,
                     decrange=decrange_default,
                     suptitle=suptitle,
                     color='red',
                     filename=filename,
                     colnames_radec=colnames_radec,
                     plotfile_suffix=plotfile_suffix)


    # key=raw_input("Enter any key to continue: ")

    # VIKING
    infile_viking = inpath + filename_viking
    viking = Table.read(infile_viking)
    viking.meta['filename'] = infile_viking
    colnames_radec = colnames_radec_viking
    if units_radec_viking == 'rad':
        viking[colnames_radec[0]] = np.rad2deg(viking[colnames_radec[0]])
        viking[colnames_radec[1]] = np.rad2deg(viking[colnames_radec[1]])

    data = viking
    print('Number of rows', len(data))
    data.info()
    data.info('stats')

    filename = filename_viking
    plotfile_suffix='viking'
    explore_data(data,
                 rarange=rarange_default,
                 decrange=decrange_default,
                 suptitle=None,
                 color='red',
                 filename=filename,
                 colnames_radec=colnames_radec,
                 plotfile_suffix=plotfile_suffix)


    # magnitude columns: ZAPERMAG3
    waveband_names_viking =['Z', 'Y', 'J', 'H', 'KS']
    magtype_viking = 'APERMAG3'
    magtype = magtype_viking
    # cycle through wavebands
    for waveband in waveband_names_viking:
        itest = (viking[waveband + magtype] > 0)
        data = viking[itest]
        plotfile_suffix = 'viking' + '_' + waveband
        suptitle = plotfile_suffix
        explore_data(data,
                     rarange=rarange_default,
                     decrange=decrange_default,
                     suptitle=suptitle,
                     color='red',
                     filename=filename,
                     colnames_radec=colnames_radec,
                     plotfile_suffix=plotfile_suffix)





    # xxl
    infile_xxl_sourcelist = inpath + filename_xxl_sourcelist
    xxl_sourcelist = Table.read(infile_xxl_sourcelist)
    xxl_sourcelist.meta['filename'] = infile_xxl_sourcelist
    print('Number of rows', len(xxl_sourcelist))
    data = xxl_sourcelist
    data.info()
    data.info('stats')

    filename = filename_xxl_sourcelist
    plotfile_suffix='xxl_sourcelist'
    explore_data(data,
                 rarange=rarange_default,
                 decrange=decrange_default,
                 suptitle=None,
                 filename=filename,
                 colnames_radec=colnames_radec_vhs,
                 plotfile_suffix=plotfile_suffix,
                 color='black')

    # xxl_redshifts needs join on UXID to get positions
    XXLRedshifts = False
    if XXLRedshifts:
        infile_xxl_redshifts = inpath + filename_xxl_redshifts
        xxl_redshifts = Table.read(infile_xxl_redshifts)
        xxl_redshift.meta['filename'] = infile_xxl_redshifts
        print('Number of rows', len(xxl_redshifts))
        data = xxl_redshifts
        data.info()
        data.info('stats')

        filename = filename_xxl_redshifts
        plotfile_suffix='xxl_redshifts'
        explore_data(data,
                     rarange=rarange_default,
                     decrange=decrange_default,
                     suptitle=None,
                     filename=filename,
                     colnames_radec=colnames_radec_vhs,
                     plotfile_suffix=plotfile_suffix,
                     color='black')


    # sys.exit()


    # xmatch GALEX and Gaia to check astrometry
    infile_galex = inpath + filename_galex
    galex = Table.read(infile_galex)
    galex.meta['filename'] = infile_galex
    print('Number of rows', len(galex))
    galex.info()

    infile_gaia = inpath + filename_gaia
    gaia = Table.read(infile_gaia)
    gaia.meta['filename'] = infile_gaia
    print('Number of rows', len(gaia))
    gaia.info()

    data1 = galex
    colnames_radec1 = colnames_radec_galex
    data2 = gaia
    colnames_radec2 = colnames_radec_gaia

    plotfile_label='galex_xnn_gaia'
    idx1, rsep = xmatch_cat(table1=data1, table2=data2,
                            colnames_radec1=colnames_radec1,
                            colnames_radec2=colnames_radec2,
                            stats=True,
                            nthneighbor=1)

    # reverse the direction
    data1 = gaia
    colnames_radec1 = colnames_radec_gaia

    data2 = galex
    colnames_radec2 = colnames_radec_galex

    rmax = 90.0
    plotfile_label='gaia_xnn_galex'
    idx1, rsep = xmatch_cat(data1=data1, data2=data2,
                              colnames_radec1=colnames_radec1,
                              colnames_radec2=colnames_radec2,
                              units_radec1=['degree', 'degree'],
                              units_radec2=['degree', 'degree'],
                              rmax=rmax,
                              plotfile_label=plotfile_label,
                              nthneighbor=1)



    # Gaia internal checks
    infile_gaia = inpath + filename_gaia
    gaia = Table.read(infile_gaia)
    gaia.meta['filename'] = infile_gaia
    print('Number of rows', len(gaia))
    gaia.info()
    gaia.info('stats')
    data = gaia
    ndata = len(data)

    colname_ra = colnames_radec_gaia[0]
    colname_dec = colnames_radec_gaia[1]

    suptitle = filename_gaia
    plotfile_label = 'gaia'
    rmax = 120.0
    markersize = min(2.0, 100000.0/ndata)
    idx_xmatch = explore_surveys_checkplots(
                     data, rmax=rmax,
                     colname_ra=colname_ra, colname_dec=colname_dec,
                     plotfile_label=plotfile_label,
                     suptitle=suptitle)

    plotfile = 'plot_selfxmatch_gaia_zoom.png'
    rmax = 10.0
    idx_xmatch = xmatch_selfcheck(
                     data=data, rmax=rmax,
                     suptitle=suptitle,
                     colnames_radec=[colname_ra, colname_dec],
                     units_radec=['degree', 'degree'],
                     plotfile=plotfile)


    # GALEX internal checks
    infile_galex = inpath + filename_galex
    galex = Table.read(infile_galex)
    galex.meta['filename'] = infile_galex
    data = galex
    ndata = len(data)
    print('Number of rows', len(data))
    data.info()
    data.info('stats')

    colname_ra = 'RAJ2000'
    colname_dec = 'DEJ2000'
    suptitle = filename_galex
    plotfile_label = 'galex'

    ra = data[colname_ra]
    dec = data[colname_dec]

    rarange= [np.min(ra), np.max(ra)]
    decrange = [np.min(dec), np.max(dec)]
    suptitle = filename_galex
    plotfile = 'plot_radec_galex.png'
    ndata = len(ra)
    markersize = min(1.0, 100000.0/ndata)
    plot_radec(ra, dec, plotfile=plotfile,
               markersize=markersize,
               rarange=rarange, decrange=decrange,
               suptitle=suptitle)

    plt.close()


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

    data2 = gaia
    colnames_radec2 = colnames_radec_gaia

    plotfile_label='xxl_source_xnn_gaia'
    idx1, rsep = xmatch_cat(data1=data1, data2=data2,
                              colnames_radec1=colnames_radec1,
                              colnames_radec2=colnames_radec2,
                              xmatch_rmax=None,
                              plot_rmax=None,
                              separations=False,
                              plotfile_label=plotfile_label,
                              nthneighbor=1)

    key=raw_input("Enter any key to continue: ")



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
    xmatch_checkplot(ra1, dec1, ra2, dec2,
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
    xmatch_selfcheck(data=data, rmax=rmax, suptitle=suptitle,
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
    xmatch_selfcheck(data=data, rmax=rmax, suptitle=suptitle,
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
    xmatch_selfcheck(data=data, rmax=rmax, suptitle=suptitle,
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
