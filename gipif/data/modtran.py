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

import os
import datetime
import commands
from gipif.utils import List2File


def atmospheric_model(doy, lat):
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
        else:
            model = 5
    else:
        if summer:
            model = 2
        else:
            model = 3
    return model


class MODTRAN():
    # hard-coded options
    filterfile = True
    _workdir = 'modtran'
    _datadir = '/usr/local/modtran/DATA'

    def __init__(self, meta, date_time, lat, lon, merraprofile=False):
        self.meta = meta

        self.lat = lat
        self.lon = lon
        self.datetime = date_time
        seconds = (date_time.second + date_time.microsecond/1000000.)/3600.
        self.dtime = self.datetime.hour + self.datetime.minute/60.0 + seconds
        self.julianday = (date_time - datetime.datetime(date_time.year, 1, 1)).days + 1

        self.model = atmospheric_model(self.julianday, lat)

        #fout = open('atm.txt','w')
        #fout.write('{:>5}{:>20}{:>20}\n'.format('Band','%T','Radiance'))

        # TODO - Update to use tmp directory
        if not os.path.exists(self._workdir):
            os.makedirs(self._workdir)
        pwd = os.getcwd()
        os.chdir(self._workdir)

        # Create link to MODTRAN data dir
        if not os.path.lexists('DATA'):
            os.symlink(self._datadir, 'DATA')

        if merraprofile:
            mprofile = atmprofile(lat, lon, date)
            height = mprofile['height']
            temp = mprofile['temp']
            humidity = mprofile['humidity']
            ozone = mprofile['ozone']
            self.atmprofile = []
            for i in reversed(range(0, len(height))):
                self.atmprofile.append(self.card2c1(P=height[i], T=temp[i], H2O=humidity[i]*1000, O3=ozone[i]*1000))
            #from pprint import pprint
            #pprint(self.atmprofile)
        else:
            self.atmprofile = None

        # Generate MODTRAN input files

        # Determine if radiance or transmittance mode
        rootnames = self.addband()
        List2File(rootnames, 'mod5root.in')

        # run output and get results
        modout = commands.getstatusoutput('modtran')
        self.output = self.readoutput()

        # Change back to original directory
        os.chdir(pwd)

    def readoutput(self):
        bandnum = self.meta['bandnum']
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
            return [trans, Lu, Ld]
        except IOError as e:
            #print 'No MODTRAN data for band ',bandnum
            return

    def addband(self):
        bandnum = self.meta['bandnum']
        rootname1 = 'band' + str(bandnum)
        wvlen1 = self.meta['wvlen1']
        wvlen2 = self.meta['wvlen2']
        if self.meta['wvlen'] < 3:
            """ Run in transmittance mode for visible bands """
            mode = 4
            fwhm = 0.001
        else:
            """ Run in radiance mode for MWIR and LWIR bands """
            mode = 2
            fwhm = 0.1
        # Write tape5 file
        self.tape5(rootname1, mode, wvlen1, wvlen2, fwhm)
        if mode == 2:
            rootname2 = rootname1 + 'Ld'
            self.tape5(rootname2, mode, wvlen1, wvlen2, fwhm, surref=1, h1=0.001)
            return (rootname1, rootname2)
        else:
            return (rootname1,)

    def tape5(self, fname, mode, wvlen1, wvlen2, fwhm, surref=0, h1=100):
        f = open(fname+'.tp5', 'w')
        f.write(self.card1(mode=mode, surref=surref)+'\n')
        f.write(self.card1a()+'\n')
        if self.filterfile:
            f.write(self.card1a3()+'\n')
        f.write(self.card2()+'\n')
        if self.atmprofile is not None:
            f.write(self.card2c(len(self.atmprofile))+'\n')
            for i in self.atmprofile:
                f.write(i+'\n')
        f.write(self.card3(h1=h1)+'\n')
        f.write(self.card3a1()+'\n')
        f.write(self.card3a2()+'\n')
        f.write(self.card4(wvlen1, wvlen2, fwhm)+'\n')
        f.write(self.card5()+'\n')
        f.close()

    def card1(self, mode, surref):
        MODTRN = 'M'  # 'C' for correlated k
        card = ('{MODTRN:1}{SPEED:1}{BINARY:1}{LYMOLC:1}{MODEL:1d}{T_BEST:1}{ITYPE:4d}{IEMSCT:5d}{IMULT:5d}'
                '{M1:5d}{M2:5d}{M3:5d}{M4:5d}{M5:5d}{M6:5d}{MDEF:5d}{I_RD2C:5d} {NOPRNT:4d}{TPTEMP:8.4f}{SURREF:>7}')
        sm = self.model
        if self.atmprofile is not None:
            model = 8
            rd2c = 1
        else:
            model = sm
            rd2c = 0
        return card.format(MODTRN=MODTRN,SPEED='',BINARY='',LYMOLC='',MODEL=model,T_BEST='',ITYPE=3,IEMSCT=mode,IMULT=-1,
                    M1=sm,M2=sm,M3=sm,M4=sm,M5=sm,M6=sm,MDEF=0,I_RD2C=rd2c,NOPRNT=-1,TPTEMP=0.001,SURREF=surref)

    def card1a(self):
        if self.filterfile:
            LFLTNM = 'T'
        else:
            LFLTNM = 'F'
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
        return card.format(IPARM=1,IPH=2,IDAY=self.julianday,ISOURC=1)

    def card3a2(self):
        card = ('{PARM1:10.3f}{PARM2:10.3f}{PARM3:10.3f}{PARM4:10.3f}{TIME:10.3f}{PSIPO:10.3f}{ANGLEM:10.3f}{G:10.3f}')
        return card.format(PARM1=self.lat,PARM2=self.lon,PARM3=0,PARM4=0,TIME=self.dtime,PSIPO=0,ANGLEM=0,G=0)

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

