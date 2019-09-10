""" 
Given a specific position and time and search radius
and search window (in time),
check IPN, Fermi, and Swift to see whether any GRBs
were detected.

Note that this assumes that the year is 20XX.
"""

import numpy as np
import requests
import lxml.html as lh
import pandas as pd
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

    # Make sure that each day is represented
    window_grid = Time(
            np.arange(window[0].jd, window[-1].jd+1, 0.5), format='jd')
    searchstr = [get_searchstr(t) for t in window_grid]
    searchstr = np.unique(np.array(searchstr))

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
        print(burst_datetime)
        if np.logical_and(
                burst_datetime >= window[0], burst_datetime <= window[-1]):
            final_set.append(l)

    print("There are %s bursts in the %s-day window" %(
        len(final_set), window[-1]-window[0]))

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

    the html scraper borrows heavily from
    https://towardsdatascience.com/web-scraping-html-tables-with-python-c9baba21059
    """
    print('\n')
    print("CONDUCTING SEARCH OF FERMI SUBTHRESHOLD CATALOG")
    url = "https://gcn.gsfc.nasa.gov/fermi_gbm_subthresh_archive.html"
    page = requests.get(url)
    doc = lh.fromstring(page.content)
    tr_elements = doc.xpath('//tr')

    # the first element is notes
    # the second element is the header
    col = []
    for ii,t in enumerate(tr_elements[1]):
        name = t.text_content()
        col.append((name,[]))

    # the third element onwards are rows of the table
    for j in range(2, len(tr_elements)):
        T = tr_elements[j]
        i = 0
        for t in T.iterchildren():
            data = t.text_content()
            col[i][1].append(data)
            i+=1

    # create dataframe
    Dict={title:column for (title,column) in col}
    df=pd.DataFrame(Dict)


    # Now, go through the table and check whether each burst
    # is consistent with the time interval.
    # If so, calculate the RA and Dec offset
    for index,row in df.iterrows():
        date = '20' + str('-'.join(row['Date'].split('/')))
        time = row['Time UT']
        datetime = Time('%sT%s' %(date,time), format='isot')
        if np.logical_and(datetime > window[0], datetime < window[-1]):
            if row['Rel'] != 2:
                ra = row['RA(J2000)[deg]']
                dec = row['Dec(J2000)[deg]']
                epos = row['Error[deg]']
                print("Subthreshold burst at %s, %s with error %s" %
                        (ra, dec, epos))
                c2  = SkyCoord(ra, dec, unit='deg')
                dist = c.separation(c2).degree
                print("Distance is %s deg from source" %dist)
                return ""
    print("No Fermi subthreshold bursts found")

######################################

def swift():
    """ 
    check the Swift GRB table
    https://swift.gsfc.nasa.gov/archive/grb_table/

    I downloaded both the 2018 and 2019 tables,
    assuming that I won't need to search outside that window.
    """
    print("\n")
    print("CONDUCTING SEARCH OF SWIFT/BAT CATALOG")
    dat = np.loadtxt("swift_grb_2018.txt", dtype=str, skiprows=1) 
    nrows = dat.shape[0]
    year = ["20%s-%s-%s" %(i[0:2],i[2:4],i[4:6]) for i in dat[:,0]]
    time = dat[:,1]
    ra = dat[:,2].astype(float)
    dec = dat[:,3].astype(float)
    epos = dat[:,4].astype(float)
    date = [Time(
        '%sT%s' %(year[i],time[i]), format='isot') for i in np.arange(nrows)]

    dat = np.loadtxt("swift_grb_2019.txt", dtype=str, skiprows=1) 
    nrows = dat.shape[0]
    year = ["20%s-%s-%s" %(i[0:2],i[2:4],i[4:6]) for i in dat[:,0]]
    time = dat[:,1]
    date2 = [Time(
        '%sT%s' %(year[i],time[i]), format='isot') for i in np.arange(nrows)]

    epos_temp = dat[:,4]
    epos_temp[epos_temp=='n/a'] = '0'

    ra = np.hstack((ra, dat[:,2].astype(float)))
    dec = np.hstack((dec, dat[:,3].astype(float)))
    epos = np.hstack((epos, epos_temp.astype(float)))
    date.extend(date2)
    date =np.array(date)
    
    keep = np.logical_and(date >= window[0], date <= window[-1])
    if sum(keep) == 0:
        print("No Swift/BAT bursts during this window")
    else:
        print("%s Swift/BAT bursts during this window" %sum(keep))
        for i in np.where(keep)[0]:
            print("Position is %s, %s with error %s deg" %(
                ra[i],dec[i],epos[i]))
            c2 = SkyCoord(ra[i], dec[i], unit='deg')
            print("Separation is %s deg" %c.separation(c2).degree)


if __name__=="__main__":
    # Koala search parameters, for testing
    ra = 29.880046
    dec = 23.845460
    c = SkyCoord(ra, dec, unit='deg')
    end_time = Time(2458490.6227, format="jd")
    start_time = Time(2458488.6602, format="jd")
    window = [start_time, end_time]

    ipn()
    fermi()
    fermi_subthreshold()
    swift()
