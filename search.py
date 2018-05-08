""" 
This script takes in an RA, Dec, and Time.
It searches the GBM history to see whether the detector was pointing that way
at the time. 
"""

import numpy as np
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
    # However, in practice, you only need the year, month, and day
    # as a searchstring YYMMDD
    gtiflag = GBMgeo.checkGTI(met)


if __name__=="__main__":
    ra = 221.491713
    dec = 14.993165
    # burst time of iPTF14yb
    t = Time('2014-02-26T10:02:57', format='isot', scale='utc')

