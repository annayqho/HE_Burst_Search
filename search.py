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


if __name__=="__main__":
    ra = 221.491713
    dec = 14.993165
    # burst time of iPTF14yb
    t = Time('2014-02-26T10:02:57', format='isot', scale='utc')
    search_gbm(ra, dec, t)
