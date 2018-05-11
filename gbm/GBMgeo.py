# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 14:51:44 2016

@author: burnse
"""


from astropy.io import fits
from astropy.time import Time
import astropy.units as u
from os import listdir
from os.path import isfile, join
import numpy as np
import math
#import healpy as hp
#import JointLoc

#import math
import urllib.request
import datetime
#from astropy.io import fits
from datetime import timedelta
from urllib.request import HTTPError
#import getGBMdata
import glob
from ftplib import FTP
import os





poshistpath = '/Users/annaho/Dropbox/Projects/Research/HE_Burst_Search/'

naitable = """
 n0  45.9  20.6
 n1  45.1  45.3
 n2  58.4  90.2
 n3 314.9  45.2
 n4 303.2  90.3
 n5   3.4  89.8
 n6 224.9  20.4
 n7 224.6  46.2
 n8 236.6  90.0
 n9 135.2  45.6
 na 123.7  90.4
 nb 183.7  90.3
"""
naitable = [lined.strip().split() for lined in naitable.strip().split('\n')]
bgotable = [['b0', 0, 90], ['b1', 180, 90]]
dettable = naitable+bgotable


def getAngOff(ara, adec, bra, bdec):
    ara = ara*math.pi/180.
    adec = adec*math.pi/180.
    bra = bra*math.pi/180.
    bdec = bdec*math.pi/180.
    tempterm = math.sin(adec)*math.sin(bdec)+math.cos(adec)*math.cos(bdec)*math.cos(bra-ara)
    if tempterm > 1.0 and tempterm < 1.0000000000000005:
        tempterm = 1.0
    try:
        AngOff = math.acos(tempterm)*180./math.pi
        return AngOff
    except ValueError:
        print('Problem with the offset angle calculation')
 
def getData(MET, putdir=None, basedir='', getTTEflag=False, getPOSflag=False, getCSPECflag=False, getCTIMEflag=False, getAllflag=False):

    if getAllflag == True:
        getTTEflag = True
        getPOSflag = True
        getCSPECflag = True
        getCTIMEflag = True
        
    utc = fermi2utc(MET)
    yymmddpath = utc.strftime('%y')+'/'+utc.strftime('%m')+'/'+utc.strftime('%d')+'/'
    
    if putdir == None:     
        putdir = basedir+'/20'+yymmddpath+'current/'
    if not os.path.exists(putdir):
        os.makedirs(putdir)
       
    ftp = FTP('legacy.gsfc.nasa.gov')
    ftp.login()
    ftpdir = 'fermi/data/gbm/daily/20'+yymmddpath+'current/'
    ftp.cwd(ftpdir)
    foldfiles = ftp.nlst()
    
    try:
        if getPOSflag == True:
            tempfilename = [f for f in foldfiles if 'poshist' in f][0]
            if len(tempfilename) > 0:
                filename = [f for f in foldfiles if 'poshist' in f][0]
                putfold = putdir+filename
                if not os.path.isfile(putfold):
                    putfile = open(putfold, 'wb')
                    ftp.retrbinary("RETR " + filename, putfile.write)
                    putfile.close()
                    
        if getCTIMEflag == True:
            flist = [f for f in foldfiles if 'ctime' in f]
            for fn in flist:
                putfold = putdir+fn
                if not os.path.isfile(putfold):
                    putfile = open(putfold, 'wb')
                    ftp.retrbinary("RETR " + fn, putfile.write)
                    putfile.close()
        
        if getCSPECflag == True:
            flist = [f for f in foldfiles if 'cspec' in f]
            for fn in flist:
                putfold = putdir+fn
                if not os.path.isfile(putfold):
                    putfile = open(putfold, 'wb')
                    ftp.retrbinary("RETR " + fn, putfile.write)
                    putfile.close()

        if getTTEflag == True:
            
            dayfrac = int(((utc.hour*60.0+utc.minute)*60.0+utc.second)/86400.*1000.)
            
            flist = [f for f in foldfiles if 'tte' in f]
            
            tlist = []
            for thing in flist:
                pt = float(thing[17:20])
                if pt not in tlist:
                    tlist.append(pt)
            
            comp = 1000
            for t in tlist:
                if dayfrac > t:
                    if dayfrac - t < comp:
                        comp = dayfrac-t
            goodfrac = str(int(dayfrac-comp))
            
            goodstr = utc.strftime('%y%m%d')+str(goodfrac).zfill(3)
            flist = [x for x in flist if goodstr in x]            
                
            for fn in flist:
                putfold = putdir+fn
                if not os.path.isfile(putfold):
                    putfile = open(putfold, 'wb')
                    ftp.retrbinary("RETR " + fn, putfile.write)
                    putfile.close()

        ftp.quit()    
    except IndexError:
        print('Index Error in getGBMdata')
        ftp.quit()

def getDetAngles(MET, inra, indec, poshistloc=poshistpath):
    sdt = datetime.datetime(2001,1,1,0,0,0)+timedelta(seconds=MET)
    
    calendate = sdt.strftime('%y')+sdt.strftime('%m')+sdt.strftime('%d')
    
    poshists = glob.glob(poshistloc+'glg_poshist_all_'+calendate+'_*.fit')
    if len(poshists) == 0:
        poshistflag = False
        i = 0
        while poshistflag == False and i < 5:
            fname = 'glg_poshist_all_'+calendate+'_v0'+str(i)+'.fit'
            storespot = poshistloc+fname
            startstring = 'http://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/'
            datestring = str(sdt.year)+'/'+sdt.strftime('%m')+'/'+sdt.strftime('%d')+'/current/'
            url = startstring+datestring+fname
            try:
                urllib.request.urlretrieve(url, storespot)
                poshistflag = True
            except HTTPError:
                i = i+1
            except FileNotFoundError:
                poshistloc = 'D:/OneDrive/Work/V404/background/poshistfiles/'
                storespot = poshistloc+fname
                try:
                    urllib.request.urlretrieve(url, storespot)
                    poshistflag = True
                except HTTPError:
                    i = i+1

    else:
        storespot = poshists[0]
        poshistflag = True
    
    if poshistflag == False:
        print('Couldn\'t find a poshist file in '+url)
    else:
        with fits.open(storespot) as hdulist:
            totdata = hdulist[1].data
        timecol = totdata.field(0)
        
        rownum = np.abs(timecol-MET).argmin()
        
        ####Gets variables necessary for calculating Earth offset angle (to determine if occulted)  
        sc_quat0 = totdata.field('QSJ_1')[rownum]
        sc_quat1 = totdata.field('QSJ_2')[rownum]
        sc_quat2 = totdata.field('QSJ_3')[rownum]
        sc_quat3 = totdata.field('QSJ_4')[rownum]
        
        ####Calculates x, y, z position from quaternions
        scx0 = sc_quat0*sc_quat0-sc_quat1*sc_quat1-sc_quat2*sc_quat2+sc_quat3*sc_quat3
        scx1 = 2*(sc_quat0*sc_quat1 + sc_quat3*sc_quat2)
        scx2 = 2*(sc_quat0*sc_quat2 - sc_quat3*sc_quat1)
        scy0 = 2*(sc_quat0*sc_quat1 - sc_quat3*sc_quat2)
        scy1 = -sc_quat0*sc_quat0+sc_quat1*sc_quat1-sc_quat2*sc_quat2+sc_quat3*sc_quat3
        scy2 = 2*(sc_quat1*sc_quat2 + sc_quat3*sc_quat0)
        scz0 = 2*(sc_quat0*sc_quat2 + sc_quat3*sc_quat1)
        scz1 = 2*(sc_quat1*sc_quat2 - sc_quat3*sc_quat0)
        scz2 = -sc_quat0*sc_quat0-sc_quat1*sc_quat1+sc_quat2*sc_quat2+sc_quat3*sc_quat3
        
        dtorad = 180./math.pi
        fra = inra/dtorad
        fdec = indec/dtorad
        
        source_pos0 = math.cos(fdec)*math.cos(fra)
        source_pos1 = math.cos(fdec)*math.sin(fra)
        source_pos2 = math.sin(fdec)
        
        sc_source_pos0 = scx0*source_pos0+scx1*source_pos1+scx2*source_pos2
        sc_source_pos1 = scy0*source_pos0+scy1*source_pos1+scy2*source_pos2
        sc_source_pos2 = scz0*source_pos0+scz1*source_pos1+scz2*source_pos2
        
        scazi = math.atan2(sc_source_pos1, sc_source_pos0)*180.0/math.pi
        if scazi < 0:
            scazi = 360+scazi
        sczen = math.acos(sc_source_pos2/1.0)*180.0/math.pi
        
        outtabletemp = []
        for thing in naitable:
            thra = float(thing[1])
            thdec = 90-float(thing[2])
            angoff = getAngOff(thra, thdec, scazi, 90-sczen)
            outtabletemp.append([thing[0], angoff])
        outtable = []
                
        bgotable = [['b0', 0, 0], ['b1', 180, 0]]
        bgo0angoff = getAngOff(float(bgotable[0][1]), float(bgotable[0][2]), scazi, 90-sczen)
        bgo1angoff = getAngOff(float(bgotable[1][1]), float(bgotable[1][2]), scazi, 90-sczen)
        outtable.append([bgotable[0][0], bgo0angoff])
        outtable.append([bgotable[1][0], bgo1angoff])
        
        outfin = outtabletemp+outtable
        return outfin









def getPointing(MET, poshistloc=poshistpath):
    
    sdt = datetime.datetime(2001,1,1,0,0,0)+timedelta(seconds=MET)
    
    calendate = sdt.strftime('%y')+sdt.strftime('%m')+sdt.strftime('%d')
    
    poshists = glob.glob(poshistloc+'glg_poshist_all_'+calendate+'_*.fit')
    if len(poshists) == 0:
        poshistflag = False
        i = 0
        while poshistflag == False and i < 5:
            fname = 'glg_poshist_all_'+calendate+'_v0'+str(i)+'.fit'
            storespot = poshistloc+fname
            startstring = 'http://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/'
            datestring = str(sdt.year)+'/'+sdt.strftime('%m')+'/'+sdt.strftime('%d')+'/current/'
            url = startstring+datestring+fname
            try:
                urllib.request.urlretrieve(url, storespot)
                poshistflag = True
            except HTTPError:
                i = i+1
    else:
        storespot = poshists[0]
        poshistflag = True
        
    
    if poshistflag == False:
        print('Couldn\'t find a poshist file in '+url)
    else:
        with fits.open(storespot) as hdulist:
            totdata = hdulist[1].data
        timecol = totdata.field(0)
        
        rownum = np.abs(timecol-MET).argmin()
        
        ####Gets variables necessary for calculating Earth offset angle (to determine if occulted)  
        sc_quat0 = totdata.field('QSJ_1')[rownum]
        sc_quat1 = totdata.field('QSJ_2')[rownum]
        sc_quat2 = totdata.field('QSJ_3')[rownum]
        sc_quat3 = totdata.field('QSJ_4')[rownum]
        
        ####Calculates x, y, z position from quaternions
        scx0 = sc_quat0*sc_quat0-sc_quat1*sc_quat1-sc_quat2*sc_quat2+sc_quat3*sc_quat3
        scx1 = 2*(sc_quat0*sc_quat1 + sc_quat3*sc_quat2)
        scx2 = 2*(sc_quat0*sc_quat2 - sc_quat3*sc_quat1)
        scy0 = 2*(sc_quat0*sc_quat1 - sc_quat3*sc_quat2)
        scy1 = -sc_quat0*sc_quat0+sc_quat1*sc_quat1-sc_quat2*sc_quat2+sc_quat3*sc_quat3
        scy2 = 2*(sc_quat1*sc_quat2 + sc_quat3*sc_quat0)
        scz0 = 2*(sc_quat0*sc_quat2 + sc_quat3*sc_quat1)
        scz1 = 2*(sc_quat1*sc_quat2 - sc_quat3*sc_quat0)
        scz2 = -sc_quat0*sc_quat0-sc_quat1*sc_quat1+sc_quat2*sc_quat2+sc_quat3*sc_quat3
        
        latboresighttheta = 0.0
        latboresightphi = 0.0
        sc_source_pos0 = math.sin(latboresighttheta)*math.cos(latboresightphi)
        sc_source_pos1 = math.sin(latboresighttheta)*math.sin(latboresightphi)
        sc_source_pos2 = math.cos(latboresighttheta)
    
        source_x = scx0*sc_source_pos0+scy0*sc_source_pos1+scz0*sc_source_pos2
        source_y = scx1*sc_source_pos0+scy1*sc_source_pos1+scz1*sc_source_pos2
        source_z = scx2*sc_source_pos0+scy2*sc_source_pos1+scz2*sc_source_pos2
        
        
        source_theta = math.acos(source_z/math.sqrt(source_x*source_x+source_y*source_y+source_z*source_z))
        source_phi = math.atan2(source_y, source_x)
        
        lat_ra = source_phi*180.0/math.pi
        if lat_ra < 0.0:
            lat_ra = lat_ra+360.0
        lat_dec = 90.0-source_theta*180.0/math.pi
        
        
        return lat_ra, lat_dec





def metToYYMMDD(MET):
    """ Convert MET (Mission Elapsed Time)
    to the searchstring needed for checkGTI 
    
    Parameters
    ----------
    MET: Mission Elapsed Time
    """
    # cMET: Mission Elapsed Time
    # defined as the time since this reference below
    t0 = Time('2001-01-01T00:00:00')
    dateobs = str((t0+(cMET)*u.s).value)[0:19]
    diryear = str(dateobs[0:4])
    dirmonth = str(dateobs[5:7])
    dirday = str(dateobs[8:10])
    searchstring = diryear[2:]+dirmonth+dirday
    return searchstring




def checkGTI(yymmdd, poshistbase=poshistpath):
    """ Check if mission was operating at the time

    Parameters
    ----------
    yymmdd: string of date to check
    """
    
    if not os.path.isdir(poshistbase):
        poshistbase = 'D:/OneDrive/Work/V404/background/poshistfiles/'

    phf = [ f for f in listdir(poshistbase) if isfile(join(poshistbase,f)) and yymmdd in f]
    
    if not os.path.isfile(poshistbase+phf[0]):    
        getData(cmet, putdir=poshistbase, getPOSflag=True)

    with fits.open(poshistbase+phf[0]) as hdulist:
        trigtotdata = hdulist[1].data
    times = trigtotdata.field("SCLK_UTC")
    tidx = (np.abs(times-cMET)).argmin()
    
    saaflag = trigtotdata[tidx][-1]
    if saaflag == 1.0:
        return True
    else:
        return False




def getEarthCenter(cMET, poshistbase = poshistpath):
    
    
    if not os.path.isdir(poshistbase):
        poshistbase = 'D:/OneDrive/Work/V404/background/poshistfiles/'
        
    t0 = Time('2001-01-01T00:00:00')
    dateobs = str((t0+(cMET)*u.s).value)[0:19]
    diryear = str(dateobs[0:4])
    dirmonth = str(dateobs[5:7])
    dirday = str(dateobs[8:10])
    searchstring = diryear[2:]+dirmonth+dirday
    phf = [ f for f in listdir(poshistbase) if isfile(join(poshistbase,f)) and searchstring in f]
    
    with fits.open(poshistbase+phf[0]) as hdulist:
        trigtotdata = hdulist[1].data
    times = trigtotdata.field("SCLK_UTC")
    tidx = (np.abs(times-cMET)).argmin() 
    
    sc_x = -trigtotdata.field("POS_X")[tidx]
    sc_y = -trigtotdata.field("POS_Y")[tidx]
    sc_z = -trigtotdata.field("POS_Z")[tidx]
    sc_theta = math.acos(sc_z/math.sqrt(sc_x*sc_x + sc_y*sc_y + sc_z*sc_z))
    sc_phi = math.atan2(sc_y, sc_x)
    sc_ra = sc_phi*180.0/math.pi
    if sc_ra < 0:
        sc_ra = 360.0+sc_ra
    sc_dec = 90.0-sc_theta*180.0/math.pi
    
    return sc_ra, sc_dec
    
    
#def getlist(cMET, NSIDE=512):
#    sra, sdec = getEarthCenter(cMET)
#    x,y,z = JointLoc.spher2cart(sra*math.pi/180.0, (90.0-sdec)*math.pi/180.0)
#    cutlist = hp.query_disc(NSIDE, [x,y,z], 66.0*math.pi/180.0, nest=True)
#    return cutlist
