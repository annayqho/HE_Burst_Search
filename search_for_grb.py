""" 
Given a specific position and time and search radius, 
what GRBs happened in the last 2 days?
"""

import numpy as np
from lxml import etree
import urllib.request
import requests
import re
import subprocess
import os
from time import strptime
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u

def get_searchstr(t):
    """ useful function for the IPN catalog search """
    yy = str(t.datetime.year)[2:]
    mm = t.strftime('%b').upper()
    dd = str(t.datetime.day).zfill(2)
    searchstr = ' '.join([dd,mm,yy])
    return searchstr

######################################

def ipn():
    """ Check the online IPN catalog
    Note that this catalog ends on 30 June 2019

    Also assumes that bursts were in 2000 onward.
    Would need to modify for years with 19XX.
    """
    print("CONDUCTING SEARCH OF IPN CATALOG")

    searchstr = [get_searchstr(t) for t in window]

    # Pull out the relevant lines
    fpath = "http://www.ssl.berkeley.edu/ipn3/masterli.txt"
    fname = fpath.split('/')[-1]
    if os.path.exists(fname) is False:
        subprocess.call(['wget', 'http://www.ssl.berkeley.edu/ipn3/masterli.txt'])
    lines = np.array(open(fname, "r").readlines())
    header = lines[np.array([' DOY TIME ' in l for l in lines])][0]
    keep = np.array([l[7:16] in searchstr for l in lines])

    # Now, check each one to see if the time is correct
    final_set = []
    for l in lines[keep]:
        filtered_l = [i for i in l.split(" ") if i]
        burst_dd = filtered_l[0].split('.')[1]
        burst_mm = str(strptime(filtered_l[1], '%b').tm_mon).zfill(2)
        burst_yy = str(filtered_l[2]).zfill(2)
        burst_time = filtered_l[4]
        burst_datetime = Time(
                '20%s-%s-%sT%s' %(burst_yy,burst_mm,burst_dd,burst_time), 
                format='isot')
        if np.logical_and(burst_datetime > window[0], burst_datetime < window[-1]):
            final_set.append(l)

    print("There are %s bursts in the %s-day window" %(
        len(final_set), window_size))

    # Check which spacecraft observed these bursts
    for l in final_set:
        det_by = np.array(
                [header[i.start():i.end()] for i in re.finditer('YES', l)])
        print(det_by)

######################################

def fermi():
    """ Check the online Fermi burst catalog """
    print("\n")
    print("CONDUCTING SEARCH OF FERMI CATALOG")

    www = 'https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3query.pl'
    data = {}
    data['tablehead'] = 'name=heasarc_fermigbrst&description=Fermi GBM Burst Catalog&url=http://heasarc.gsfc.nasa.gov/W3Browse/fermi/fermigbrst.html&archive=Y&radius=180&mission=FERMI&priority=1&tabletype=Object'
    data['Time'] = '%s .. %s' %(window[0].mjd, window[-1].mjd)
    data['displaymode'] = 'PureTextDisplay'
    data['varon'] = 'name,ra,dec,trigger_time,error_radius'
    r = requests.post(url = www, data = data)
    out = np.array([i for i in r.text.split('\n') if i])
    header = [i for i in out[2].split('|')]
    ncands = len(out[3:])
    if ncands == 0:
        print("No GRBs in Fermi")
    else:
        print("Found %s in Fermi" %ncands)
        for ii in np.arange(ncands):
            cand = out[3+ii]
            vals = [i for i in cand.split('|')]
            print(vals)
            grbname = vals[1]
            grbra = vals[2]
            grbdec = vals[3]
            grbtime = vals[4]
            grbepos = vals[5].strip()
            print("%s with RA=%s, Dec=%s on t=%s with %s deg uncertainty"%(
                grbname,grbra,grbdec,grbtime,grbepos))
            c2 = SkyCoord('%s %s' %(grbra,grbdec), unit=(u.hourangle, u.deg))
            dist = c2.separation(c).degree
            print("The burst is %s deg away from the source" %dist)


def fermi_subthreshold():
    """ 
    check the fermi subthreshold notices
    https://gcn.gsfc.nasa.gov/fermi_gbm_subthresh_archive.html
    """
    fp = urllib.request.urlopen("https://gcn.gsfc.nasa.gov/fermi_gbm_subthresh_archive.html")
    mybytes = fp.read()
    mystr = mybytes.decode("utf8")
    fp.close()
    table = etree.HTML(mystr).find("body/table")
    rows = iter(table)
    headers = [col.text for col in next(rows)]
    for row in rows:
        values = [col.text for col in row]
        print(dict(zip(headers,values)))


if __name__=="__main__":
    window_size = 2

    # Koala search parameters, for testing
    ra = 30.063307 
    dec = 16.799255
    c = SkyCoord(ra, dec, unit='deg')
    time = Time(2458373.9075, format="jd")
    window = Time(
            np.linspace(time.jd-window_size, time.jd, window_size), format='jd')

    fermi_subthreshold()
