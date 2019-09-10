""" 
Given a specific position and time and search radius, 
when is the last time a GRB happened in this location?
"""

import subprocess
import os
from astropy.time import Time

# Search parameters, for testing
ra = 279.472820 
dec = 61.497984
erad = 1 # degree
time = Time(2458728.6798, format="jd")

######################################

# CHECK IPN

# You need the day, month, year


fpath = "http://www.ssl.berkeley.edu/ipn3/masterli.txt"
fname = fpath.split('/')[-1]
if os.path.exists(fname) is False:
    subprocess.call(['wget', 'http://www.ssl.berkeley.edu/ipn3/masterli.txt'])
for line in open(fname, "r"):
    dat = line.split(' ')
    dat = [i for i in dat if i]
    print(dat)
    print(len(dat))
