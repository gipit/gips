#!/usr/bin/python

import os, sys
import numpy as np
import commands
import datetime
import gippy

    #import Pysolar
    #el = Pysolar.GetAltitude(meta['lat'],meta['lon'],dt)
    #az = Pysolar.GetAzimuth(meta['lat'],meta['lon'],dt)
    #print 'el = ', el, 90 - meta['solarzenith']
    #print 'az = ', -180-az, meta['solarazimuth']

_merraroot = '/titan/data/merra/'
_moddatadir = '/usr/local/modtran/DATA'
_workdir = '/tmp'

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

def atmmodel(doy,lat):
    """ Determine atmospheric model
    1 - Tropical
    2 - Mid-Latitude Summer
    3 - Mid-Latitude Winter
    4 - Sub-Arctic Summer
    5 - Sub-Arctic Winter
    6 - US Standard Atmosphere
    """
    # Determine season
    if doy < 121 or doy > 274:
        if lat < 0:
            summer = True
        else:
            summer = False
    else: 
        if lat < 0:
            summer = False
        else:
            summer = True
    # Determine model
    if abs(lat) <= 15:
        model = 1
    elif abs(lat) >= 60:
        if summer:
            model = 4
        else: model = 5
    else:
        if summer:
            model = 2
        else: model = 3
    return model

def SixS(bandnum, meta):
    """ Run 6S atmospheric model for given spectral band """
    import agspy.data.aod.RetrieveOptDep as RetrieveOptDep
    from Py6S import SixS, Geometry, AeroProfile, Altitudes, Wavelength, GroundReflectance, AtmosCorr

    idoy = meta['JulianDay']
    model = atmmodel(idoy, meta['lat'])
    bandlocs = meta['bands']
    bandwidths = meta['bandwidths']
    dt = meta['datetime']
    wvlen1 = bandlocs[bandnum-1] - bandwidths[bandnum-1]/2.0
    wvlen2 = bandlocs[bandnum-1] + bandwidths[bandnum-1]/2.0

    # Get aerosols
    strdoy=str(idoy)
    if idoy < 10:
        strdoy = '0'+strdoy
    if idoy < 100:
        strdoy = '0'+strdoy
    aod,aod_source,aod_rmse = RetrieveOptDep.getaod(strdoy,str(dt.year),str(meta['lon']),str(meta['lat']))
    if aod < 0.0:
        aod = 0.17
        print 'NO AOD DATA AVAILABLE; USING 0.17 default' 

    s = SixS()
    s.atmos_profile = atmmodel(idoy,meta['lat'])
    s.geometry = Geometry.User()
    s.geometry.from_time_and_location(meta['lat'], meta['lon'], str(dt), meta['zenith'], meta['azimuth'])
    aero = AeroProfile.PredefinedType(AeroProfile.Continental)
    s.aero_profile = aero
    s.aot550 = aod
    s.altitudes = Altitudes()
    s.altitudes.set_target_sea_level()
    s.altitudes.set_sensor_satellite_level()
    s.wavelength = Wavelength(wvlen1,wvlen2)
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(GroundReflectance.GreenVegetation)
    s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromRadiance(1.0)

    #wavelengths, results = SixSHelpers.Wavelengths.run_wavelengths(s, numpy.arange(wvlen1, wvlen2, 0.001))
    s.run()
    #print 'coef = ', (s.outputs.coef_xa, s.outputs.coef_xb, s.outputs.coef_xc)
    out = s.outputs
    t = s.outputs.trans['global_gas'].upward #s.outputs.total_gaseous_transmittance
    Lu = s.outputs.atmospheric_intrinsic_radiance
    Ld = (s.outputs.direct_solar_irradiance + s.outputs.diffuse_solar_irradiance + s.outputs.environmental_irradiance)/np.pi
    #Xa = s.outputs.coef_xa
    #Xb = s.outputs.coef_xb
    #Xc = s.outputs.coef_xc
    #return {'t':t,'Lu':Lu,'Ld':Ld, 'Xa':Xa, 'Xb':Xb, 'Xc':Xc}
    return {'t':t,'Lu':Lu,'Ld':Ld}

#class MODTRAN(Atmosphere):
class MODTRAN():
    # hard-coded options
    filterfile = True
    workdir = 'modtran'

    def __init__(self,bandnum,meta,merraprofile=True):
        #workdir = os.path.join(workdir,'modtran')
        self.meta = meta

        self.model = atmmodel(self.meta['JulianDay'], self.meta['lat'])

        #fout = open('atm.txt','w')
        #fout.write('{:>5}{:>20}{:>20}\n'.format('Band','%T','Radiance'))

        # Update to use tmp directory
        if not os.path.exists(self.workdir): os.makedirs(self.workdir)
        pwd = os.getcwd()
        os.chdir(self.workdir)

        # Create link to MODTRAN data dir
        if not os.path.lexists('DATA'): os.symlink(_moddatadir, 'DATA')

        if merraprofile:
            mprofile = atmprofile(meta['lat'],meta['lon'],meta['datetime'])
            height = mprofile['height']
            temp = mprofile['temp']
            humidity = mprofile['humidity']
            ozone = mprofile['ozone']
            self.atmprofile = []
            for i in reversed(range(0,len(height))):
                self.atmprofile.append(self.card2c1(P=height[i],T=temp[i],H2O=humidity[i]*1000,O3=ozone[i]*1000)) 
            #from pprint import pprint
            #pprint(self.atmprofile)
        else: self.atmprofile = None

        # Generate MODTRAN input files
        bandloc = meta['bands'][bandnum-1]
        bandwidth = meta['bandwidths'][bandnum-1]
        # Determine if radiance or transmittance mode
        rootnames = self.addband(bandnum)

        mod5root = open('mod5root.in','w')
        for r in rootnames: mod5root.write(r+'\n')
        mod5root.close()
        modout = commands.getstatusoutput('modtran')
      
        self.output = self.readoutput(bandnum)

        # Change back to original directory
        os.chdir(pwd)

    def readoutput(self,bandnum):
        try:
            f = open('band'+str(bandnum)+'.chn')
            lines = f.readlines()
            f.close()
            data = lines[4+bandnum]
            # Get nominal band width in microns
            bandwidth = float(data[85:94]) / 1000
            # Convert from W/sr-cm2 to W/sr-m2-um 
            Lu = (float(data[59:72]) * 10000) / bandwidth
            trans = float(data[239:248])
            try:
                f = open('band'+str(bandnum)+'Ld.chn')
                lines = f.readlines()
                f.close()
                data = lines[4+bandnum]
                # Convert channel radiance to spectral radiance
                Ld = (float(data[59:72]) * 10000) / bandwidth
            except IOError as e:
                #print 'No downwelled radiance run'
                Ld = 0.0
            return {'t':trans,'Lu':Lu,'Ld':Ld}
        except IOError as e:
            #print 'No MODTRAN data for band ',bandnum
            return

    def addband(self,bandnum):
        rootname1 = 'band' + str(bandnum)
        bandloc = self.meta['bands'][bandnum-1]
        bandwidth = self.meta['bandwidths'][bandnum-1]
        wvlen1 = bandloc - bandwidth/2.0
        wvlen2 = bandloc + bandwidth/2.0 
        if bandloc < 3:
            """ Run in transmittance mode for visible bands """
            mode = 4
            fwhm = 0.001
        else:
            """ Run in radiance mode for MWIR and LWIR bands """
            mode = 2
            fwhm = 0.1
        # Write tape5 file
        self.tape5(rootname1,mode,wvlen1,wvlen2,fwhm)
        if mode == 2:
            rootname2 = rootname1 + 'Ld'
            self.tape5(rootname2,mode,wvlen1,wvlen2,fwhm,surref=1,h1=0.001)
            return (rootname1,rootname2)
        else:
            return (rootname1,)  

    def tape5(self,fname,mode,wvlen1,wvlen2,fwhm,surref=0,h1=100):
        f = open(fname+'.tp5','w') 
        f.write(self.card1(mode=mode,surref=surref)+'\n')
        f.write(self.card1a()+'\n')
        if self.filterfile: f.write(self.card1a3()+'\n')
        f.write(self.card2()+'\n')
        if self.atmprofile != None:
            f.write(self.card2c(len(self.atmprofile))+'\n')
            for i in self.atmprofile: f.write(i+'\n')
        f.write(self.card3(h1=h1)+'\n')
        f.write(self.card3a1()+'\n')
        f.write(self.card3a2()+'\n')
        f.write(self.card4(wvlen1, wvlen2, fwhm)+'\n')
        f.write(self.card5()+'\n')
        f.close()

    def card1(self,mode,surref):
        MODTRN = 'M' #'C' for correlated k
        card = ('{MODTRN:1}{SPEED:1}{BINARY:1}{LYMOLC:1}{MODEL:1d}{T_BEST:1}{ITYPE:4d}{IEMSCT:5d}{IMULT:5d}'
                '{M1:5d}{M2:5d}{M3:5d}{M4:5d}{M5:5d}{M6:5d}{MDEF:5d}{I_RD2C:5d} {NOPRNT:4d}{TPTEMP:8.4f}{SURREF:>7}')
        sm = self.model
        if self.atmprofile != None:
            model = 8
            rd2c = 1
        else:
            model = sm
            rd2c = 0
        return card.format(MODTRN=MODTRN,SPEED='',BINARY='',LYMOLC='',MODEL=model,T_BEST='',ITYPE=3,IEMSCT=mode,IMULT=-1,
                    M1=sm,M2=sm,M3=sm,M4=sm,M5=sm,M6=sm,MDEF=0,I_RD2C=rd2c,NOPRNT=-1,TPTEMP=0.001,SURREF=surref)

    def card1a(self):
        if self.filterfile:
            LFLTNM='T'
        else: LFLTNM='F'
        card = ('{DIS:1}{DISAZM:1}{DISALB:1}{NSTR:>3}{SFWHM:4.1f}{CO2MX:10.3f}{H2OSTR:>10}{O3STR:>10}'
                '{C_PROF:1}{LSUNFL:1} {LBMNAM:1} {LFLTNM:1} {H2OAER:1} {CDTDIR:1}{SOLCON:10.3f}')
                #'{CDASTM:1}{ASTMC:9.2f}{ASTMX:10.3f}{ASTMO:10.3f}{AERRH:10.3f}{NSSALB:10d}')
        return card.format(DIS='T',DISAZM='',DISALB='T',NSTR=8,SFWHM=0,CO2MX=380.0,H2OSTR='1.00',O3STR='1.00',
                    C_PROF='',LSUNFL=1,LBMNAM='f',LFLTNM=LFLTNM,H2OAER='f',CDTDIR='f',SOLCON=0.0)
                    #CDASTM='',ASTMC='',ASTMX='',ASTMO='',AERRH='',NSSALB=0)

    def card1a3(self):
        #card = ('{FILTNM:256} !CARD1A3')
        card = ('{FILTNM:56}')
        return card.format(FILTNM='DATA/landsat7.flt')

    def card2(self):
        card = ('{APLUS:>2}{IHAZE:3d}{CNOVAM:1}{ISEASN:4d}{ARUSS:>3}{IVULCN:2d}{ICSTL:5d}{ICLD:5d}'
                '{IVSA:5d}{VIS:10.5f}{WSS:10.5f}{WHH:10.5f}{RAINRT:10.5f}{GNDALT:10.5f}')
        return card.format(APLUS='',IHAZE=1,CNOVAM='',ISEASN=0,ARUSS='',IVULCN=0,ICSTL=3,ICLD=0,
                    IVSA=0,VIS=0,WSS=0,WHH=0,RAINRT=0,GNDALT=0)

    def card2c(self,numalt):
        card = ('{ML:5d}{IRD1:5d}{IRD2:5d}{HMODEL:>20}{REE:10.0f}{NMOLYC:5d}{E_MASS:10.0f}{AIRMWT:10.0f}')
        return card.format(ML=numalt,IRD1=0,IRD2=0,HMODEL='custom',REE=0,NMOLYC=0,E_MASS=0,AIRMWT=0)

    def card2c1(self,h=0,P=0,T=0.0,H2O=0,CO2=0,O3=0):
        """ Custom atmospheric layers """
        card = ('{ZM:10.3f}{P:10.3f}{T:10.3f}{WMOL1:10.3f}{WMOL2:10.3f}{WMOL3:10.3f}{JCHAR:>14}{JCHARX:1}{JCHARY:1}')
        return card.format(ZM=h,P=P,T=T,WMOL1=H2O,WMOL2=CO2,WMOL3=O3,JCHAR='ABC C         ',JCHARX=' ',JCHARY=' ')

    def card3(self,h1=100,alt=0,angle=180):
        card = ('{H1:10.3f}{H2:10.3f}{ANGLE:10.3f}{RANGE:10.3f}{BETA:10.3f}{RO:10.3f}     {LENN:>5}{PHI:10.3f}')
        return card.format(H1=h1,H2=alt,ANGLE=angle,RANGE=0,BETA=0,RO=0,LENN=0,PHI=0)

    def card3a1(self):
        card = ('{IPARM:>5}{IPH:>5}{IDAY:>5}{ISOURC:>5}')
        return card.format(IPARM=1,IPH=2,IDAY=self.meta['JulianDay'],ISOURC=1)

    def card3a2(self):
        card = ('{PARM1:10.3f}{PARM2:10.3f}{PARM3:10.3f}{PARM4:10.3f}{TIME:10.3f}{PSIPO:10.3f}{ANGLEM:10.3f}{G:10.3f}')
        return card.format(PARM1=self.meta['lat'],PARM2=self.meta['lon'],PARM3=0,PARM4=0,TIME=self.meta['DecimalTime'],PSIPO=0,ANGLEM=0,G=0)

    def card4(self,v1=0.4,v2=1.0,fwhm=0.002):
        """ Spectral parameters """
        dv = fwhm/2
        card = ('{V1:10.3f}{V2:10.3f}{DV:10.3f}{FWHM:10.3f}{YFLAG:1}{XFLAG:1}{DLIMIT:>8}{FLAGS:8}{MLFLX:3}{VRFRAC:10.3f}')
        return card.format(V1=v1,V2=v2,DV=dv,FWHM=fwhm,YFLAG='',XFLAG='',DLIMIT='',FLAGS="M       ",MLFLX='',VRFRAC=0)

    def card5(self):
        card = ('{IRPT:>5}')
        return card.format(IRPT=0)

    """ old plotting code
    if plot:
        #fig = plt.figure()
        #plt.title('Atmospheric Transmittance vs Wavelength (um)')

  
        results = table.Table()
        results.read_csv('tape7.scn',skip=11, delim=' ', ignorebad=True, hasheader=False) 
        #pdb.set_trace()

        wvlens = results.as_numpy_col(0)
        T = results.as_numpy_col(1)
        rad = results.as_numpy_col(2)
        fout.write('{:5d}{:20.10f}{:20.10f}\n'.format(bandnum,np.mean(T),np.mean(rad)))

        if plot:
            fig.add_subplot(4,2,bandnum)
            plt.title("Band %s" % bandnum)
            plt.plot(wvlens,T)
            plt.ylim([0.0,1.0])

        bandnum = bandnum+1
    if plot:
        plt.tight_layout()
        plt.show()
    fout.close()
    """    

# Eventually Atmosphere should take in a gippy.GeoImage
# and extract the parameters from the file itself
# to do that though, meta needs to be embedded in GeoImage
def atmosphere(bandnum,meta=None,merraprofile=False):
    vismodel = '6s'
    thmodel = 'modtran'
    #if bandnums == None:
    #    bandnums = range(1,len(meta['bands'])+1)
    # try reading atmcorr
    #fname = os.path.join(path,'atmcorr_'+vismodel+'_'+thmodel+'.txt')
    #if os.path.isfile(fname):
    #    lines = (line.rstrip('\n') for line in open(fname))
    #    atmcorr = {}
    #    for l in lines:
    #        row = l.split()
    #        atmcorr[int(row[0])] = {'t':float(row[1]), 'Lu':float(row[2]), 'Ld':float(row[3])}
    #    return atmcorr
    #else:
    curdir = os.getcwd()
    os.chdir(os.path.dirname(meta['metafilename']))
    #output = dict()
    #for i,bandloc in enumerate(meta['bandlocs']):
    #for b in bandnums:
    bandloc = meta['bands'][bandnum-1]
    # visible
    if bandloc < 3:
        if vismodel == '6s':
            atm = SixS(bandnum,meta)
        elif vismodel == 'modtran':
            raise Exception('modtran model not yet implemented for visible')
        else:
            raise Exception('unknown atm model %s' % model)
    # thermal
    else:
        if thmodel == 'modtran':
            atm = MODTRAN(bandnum,meta,merraprofile=merraprofile).output
        else:
            raise Exception('unknown atm model %s' % model)
    os.chdir(curdir)
    return gippy.Atmosphere(atm['t'],atm['Lu'],atm['Ld'])
