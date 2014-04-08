#!/usr/bin/env python
################################################################################
#    GIPPY: Geospatial Image Processing library for Python
#
#    Copyright (C) 2014 Matthew A Hanson
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <http://www.gnu.org/licenses/>
################################################################################

# TODO - WORK IN PROGRESS. CONVERT TO gippy Data subclass

from pdb import set_trace

def _interp(arr, x, y):
    """ Bilinear interpolation of 2x2 array. 0 <= (x,y) < 1.0 """
    b1 = arr[0,0]
    b2 = arr[1,0] - arr[0,0]
    b3 = arr[0,1] - arr[0,0]
    b4 = arr[0,0] - arr[1,0] - arr[0,1] + arr[1,1]
    return b1 + b2*x + b3*y + b4*x*y

def fetchmerra(date):
    product = 'MAI6NVANA.5.2.0'
    minutes = date.hour * 60 + date.minute
    if minutes < 180:
        date = date - datetime.timedelta(days=1)

    if date.year < 1993:
        stream = 1
    elif date.year < 2001:
        stream = 2
    else: stream = 3

    minorversion = ['0','1']

    for mv in minorversion:
        hdf = 'MERRA%s0%s.prod.assim.inst6_3d_ana_Nv.%s.hdf' % (stream, mv, date.strftime('%Y%m%d'))
        fname = os.path.join(_merraroot,product,date.strftime('%Y'),date.strftime('%m'),hdf)
        if not os.path.exists(fname):
            wfname = fname.replace(_merraroot,'ftp://goldsmr3.sci.gsfc.nasa.gov/data/s4pa/MERRA/')
            try:
                os.makedirs(os.path.dirname(fname))
            except:
                pass
            import urllib
            #start = datetime.datetime.now()
            try:
                urllib.urlretrieve(wfname, fname)
                break
            except:
                pass
        else: break
    return fname

def atmprofile(lat,lon,date):
    """ Fetch atmospheric profile from MERRA data """
    minutes = date.hour * 60 + date.minute

    # Determine best time to use
    if minutes < 180:
        # use previous day midnight
        timeindex = 0
        #date = date - datetime.timedelta(days=1)
    elif minutes < 540:
        timeindex = 1
    elif minutes < 900:
        timeindex = 2
    elif minutes < 1260:
        timeindex = 3
    else: timeindex = 0

    fname = fetchmerra(date)

    import nio
    f = nio.open_file(fname)

    #times = f.variables['TIME_EOSGRID'][:]
    heights = f.variables['Height_EOSGRID'][:]
    lats = f.variables['YDim_EOSGRID'][:]
    lons = f.variables['XDim_EOSGRID'][:]

    latinds = np.arange(len(lats))
    loninds = np.arange(len(lons))
    latind = np.interp(lat,lats,latinds)
    lonind = np.interp(lon,lons,loninds)
    lat0 = np.floor(latind).astype(int)
    lat1 = lat0+2
    lon0 = np.floor(lonind).astype(int)
    lon1 = lon0+2
    flat = latind-lat0
    flon = lonind-lon0

    # Extract variables
    PS = _interp(f.variables['PS'][timeindex,lat0:lat1,lon0:lon1], flat, flon)
    temp = []
    humidity = []
    ozone = []
    #delp = []
    for h in range(0,len(heights)):
        temp.append( _interp( f.variables['T'][timeindex,h,lat0:lat1,lon0:lon1], flat,flon)-273.15 )
        humidity.append( _interp( f.variables['QV'][timeindex,h,lat0:lat1,lon0:lon1], flat,flon) )
        ozone.append( _interp( f.variables['O3'][timeindex,h,lat0:lat1,lon0:lon1], flat,flon) )
        #delp.append( _interp( f.variables['DELP'][timeindex,h,lat0:lat1,lon0:lon1], flat,flon) )
    f.close()

    return {'height':heights,'temp':temp,'humidity':humidity,'ozone':ozone} #'delp':delp}

def _test_atmprofile():
    import matplotlib.pyplot as plt
    data = atmprofile(43,-72,'')
    h = data['height']

    plt.figure()

    plt.subplot(221)
    plt.plot(h,data['temp'])
    plt.xlabel('atm height (hPa)')
    plt.ylabel('temp (C)')

    plt.subplot(222)
    plt.plot(h,data['humidity'])
    plt.xlabel('atm height (hPa)')
    plt.ylabel('humidity (C)') 

    plt.subplot(223)
    plt.plot(h,data['ozone'])
    plt.xlabel('atm height (hPa)')
    plt.ylabel('ozone mixing ratio') 

    plt.subplot(224)
    plt.plot(h,data['delp'])
    plt.xlabel('atm height (hPa)')
    plt.ylabel('delp')

    plt.show()