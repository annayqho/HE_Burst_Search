""" 
This script takes in an RA, Dec, and Time.
It searches the GBM history to see whether the detector was pointing that way
at the time. 
"""

import numpy as np
import sys
import glob
from subprocess import call
import re
from astropy.time import Time
from gbm import GBMgeo
from gbm.clock import *
from swiftbat_python.swiftbat import swinfo



def search_gbm(ra, dec, t):
    """ Search the GBM history to see whether the detector was
    sensitive to a given position at a given time
    
    Parameters
    ----------
    ra: right ascension in degrees
    dec: declination in degrees
    t: time in astropy isot format
    """
    
    # The GBM package uses MET: Mission Elapsed Time
    cMET = utc2fermi(t.datetime)
    
    # Download the relevant poshist file
    yymmdd = re.sub('-', '', t.value[2:10])
    root = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/"
    date = "20" + yymmdd[0:2] + "/" + yymmdd[2:4] + "/" + yymmdd[4:6] 
    fname = "glg_poshist_all_%s_v00.fit" %yymmdd
    if glob.glob(fname):
        print("File already downloaded")
    else:
        get = root + date + "/current/glg_poshist_all_%s_v00.fit" %yymmdd
        print("Downloading the relevant poshist file")
        print(get)
        call(["wget", get])

    # Check for the time
    gtiflag = GBMgeo.checkGTI(cMET)
    if not gtiflag:
        print('Occurred during a bad time interval (likely SAA)')
    else:
        print('Occurred during a good time interval')
        Era, Edec = GBMgeo.getEarthCenter(cMET)
        angularoffset = GBMgeo.getAngOff(Era, Edec, ra, dec)
        if angularoffset > 68.0:
            print('The source was visible to GBM at this time')
        elif angularoffset > 66.0:
            print('The source was possibly visible to GBM at this time (between 66 and 68 deg offset)')
        else:
            print('The source was occulted at this time')


def search_bat(ra, dec, t):
    """ Search the BAT history to see whether the detector was 
    sensitive to a given position at a given time """
    # I think you want to use -v for 'verbose'
    # I think that you also want ot use the -t for 'timerange'
    # I thikn that you also want to use the -p for position
    # like 'RA,Dec'
    swinfo.main(
            ['2014-02-26T10:02:57', '-s 221.491713_14.993165'])
    #argv=[t.isot, 'verbose'])

    # sourcename can be
    # '123.4_-56.7' for RA_dec


if __name__=="__main__":
    ra = 221.491713
    dec = 14.993165
    # burst time of iPTF14yb
    print("Checking for iPTF14yb")
    t = Time('2014-02-26T10:02:57', format='isot', scale='utc')
    # search_gbm(ra, dec, t)
    # correct answer: offline due to SAA passage

    search_bat(ra, dec, t)

    # burst time of GRB 161228B, the iPTF17cw GRB
    # ra = 129.290
    # dec = 43.670
    # print("Checking for GRB 161228B")
    # t = Time('2016-12-28T13:15:40.0', format='isot', scale='utc')
    # search_gbm(ra, dec, t)

    
