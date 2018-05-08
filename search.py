""" 
This script takes in an RA, Dec, and Time.
It searches the GBM history to see whether the detector was pointing that way
at the time. 
"""

import numpy as np
import re
from astropy.time import Time
from gbm import GBMgeo


def search_gbm(ra, dec, t):
    """ Search the GBM history to see whether the detector was
    sensitive to a given position at a given time
    
    Parameters
    ----------
    ra: right ascension in degrees
    dec: declination in degrees
    t: time in astropy isot format
    """
    
    # The original GBM package uses MET: Mission Elapsed Time
    # I change this to accepting YYMMDD
    yymmdd = re.sub('-', '', t.value[2:10])
    # Currently can't run this part because I don't have these poshistfiles
    # gtiflag = GBMgeo.checkGTI(yymmdd)


if __name__=="__main__":
    ra = 221.491713
    dec = 14.993165
    # burst time of iPTF14yb
    t = Time('2014-02-26T10:02:57', format='isot', scale='utc')

