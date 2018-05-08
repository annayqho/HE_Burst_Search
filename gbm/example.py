# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 12:47:50 2016

@author: burnse
"""

import GBMgeo

         #[met, source ra, source dec]
data = [[463917049.391, 125.0, -4.0],
        [463927049.391, 225.0, -4.0],
        [463947049.391, 125.0, -4.0]]

for row in data:
    
    met = row[0]
    sra = row[1]
    sdec = row[2]
    
    print(met)
    
    gtiflag = GBMgeo.checkGTI(met)
    if not gtiflag:
        print('Occurred during a bad time interval (likely SAA)')
    else:
        print('Occurred during a good time interval')
        Era, Edec = GBMgeo.getEarthCenter(met)
        angularoffset = GBMgeo.getAngOff(Era, Edec, sra, sdec)
        if angularoffset > 68.0:
            print('The source was visible to GBM at this time')
        elif angularoffset > 66.0:
            print('The source was possibly visible to GBM at this time (between 66 and 68 deg offset)')
        else:
            print('The source was occulted at this time')

    print('')

