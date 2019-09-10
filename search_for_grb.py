""" 
Given a specific position and time and search radius, 
what GRBs happened in the last 2 days?
"""

import numpy as np
import requests
import re
import subprocess
import os
from astropy.time import Time

window_size = 2

# Koala search parameters, for testing
ra = 30.063307 
dec = 16.799255
erad = 1 # degree
time = Time(2458373.9075, format="jd")
window = Time(
        np.linspace(time.jd-window_size, time.jd, window_size), format='jd')

######################################

# CHECK IPN
# Note: catalog ends 30 June 2019

def get_searchstr(t):
    yy = str(t.datetime.year)[2:]
    mm = t.strftime('%b').upper()
    dd = str(t.datetime.day).zfill(2)
    searchstr = ' '.join([dd,mm,yy])
    return searchstr

searchstr = [get_searchstr(t) for t in window]

# Pull out the relevant lines
fpath = "http://www.ssl.berkeley.edu/ipn3/masterli.txt"
fname = fpath.split('/')[-1]
if os.path.exists(fname) is False:
    subprocess.call(['wget', 'http://www.ssl.berkeley.edu/ipn3/masterli.txt'])
lines = np.array(open(fname, "r").readlines())
header = lines[np.array([' DOY TIME ' in l for l in lines])][0]
keep = np.array([l[7:16] in searchstr for l in lines])
print("There are %s bursts in the last %s days" %(sum(keep), window_size))

# Check which spacecraft observed these bursts
for l in lines[keep]:
    det_by = np.array(
            [header[i.start():i.end()] for i in re.finditer('YES', l)])
    print(det_by)

######################################

# CHECK FERMI
www = 'https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3query.pl'
data = {}
data['tablehead'] = 'name=heasarc_fermigbrst&description=Fermi GBM Burst Catalog&url=http://heasarc.gsfc.nasa.gov/W3Browse/fermi/fermigbrst.html&archive=Y&radius=180&mission=FERMI&priority=1&tabletype=Object'
data['Time'] = '%s .. %s' %(window[0].mjd, window[-1].mjd)
data['displaymode'] = 'PureTextDisplay'
r = requests.post(url = www, data = data)
out = np.array([i for i in r.text.split('\n') if i])
header = [i for i in out[2].split('|')]
ncands = len(out[3:])
if ncands == 0:
    print("No GRBs in Fermi")
else:
    print("Found %s in Fermi" %ncands)
    for ii in np.arange(ncands):
        cands = out[ii]
        vals = [i for i in out[3+ii].split('|')]
        grbname = vals[1]
        grbra = vals[2]
        grbdec = vals[3]
        grbtime = vals[4]
        print("%s with RA=%s, Dec=%s on t=%s" %(grbname,grbra,grbdec,grbtime))

# subthreshold notices: 
# https://gcn.gsfc.nasa.gov/fermi_gbm_subthresh_archive.html
