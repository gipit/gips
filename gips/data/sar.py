#!/usr/bin/env python
################################################################################
#    GIPS: Geospatial Image Processing System
#
#    AUTHOR: Matthew Hanson
#    EMAIL:  matt.a.hanson@gmail.com
#
#    Copyright (C) 2014 Applied Geosolutions
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
import tarfile
import copy
import numpy

import gippy
from gips.data.core import Repository, Asset, Data
from gips.utils import File2List, List2File, RemoveFiles


class SARRepository(Repository):
    name = 'SAR'
    description = 'Synthetic Aperture Radar PALSAR and JERS-1'

    @classmethod
    def feature2tile(cls, feature):
        """ Get tile designation from a geospatial feature (i.e. a row) """
        fldindex_lat = feature.GetFieldIndex("lat")
        fldindex_lon = feature.GetFieldIndex("lon")
        lat = int(feature.GetField(fldindex_lat) + 0.5)
        lon = int(feature.GetField(fldindex_lon) - 0.5)
        if lat < 0:
            lat_h = 'S'
        else:
            lat_h = 'N'
        if lon < 0:
            lon_h = 'W'
        else:
            lon_h = 'E'
        tile = lat_h + str(abs(lat)).zfill(2) + lon_h + str(abs(lon)).zfill(3)
        return tile


class SARAsset(Asset):
    """ Single original file """
    Repository = SARRepository

    _sensors = {
        'AFBS': {'description': 'PALSAR FineBeam Single Polarization'},
        'AFBD': {'description': 'PALSAR FineBeam Dual Polarization'},
        'AWB1': {'description': 'PALSAR WideBeam (ScanSAR Short Mode)'},
        'JFBS': {'description': 'JERS-1 FineBeam Single Polarization'}
    }
    _assets = {
        '': {
            'pattern': 'KC_*.tar.gz'
        }
    }

    _defaultresolution = [0.000834028356964, 0.000834028356964]

    # launch dates for PALSAR (A) and JERS-1 (J)
    _launchdate = {'A': datetime.date(2006, 1, 24), 'J': datetime.date(1992, 2, 11)}

    _cycledates = {
        7: '20-Oct-06',
        8: '05-Dec-06',
        9: '20-Jan-07',
        10: '07-Mar-07',
        11: '22-Apr-07',
        12: '07-Jun-07',
        13: '23-Jul-07',
        14: '07-Sep-07',
        15: '23-Oct-07',
        16: '08-Dec-07',
        17: '23-Jan-08',
        18: '09-Mar-08',
        19: '24-Apr-08',
        20: '09-Jun-08',
        21: '25-Jul-08',
        22: '09-Sep-08',
        23: '25-Oct-08',
        24: '10-Dec-08',
        25: '25-Jan-09',
        26: '12-Mar-09',
        27: '27-Apr-09',
        28: '12-Jun-09',
        29: '28-Jul-09',
        30: '12-Sep-09',
        31: '28-Oct-09',
        32: '13-Dec-09',
        33: '28-Jan-10',
        34: '15-Mar-10',
        35: '30-Apr-10',
        36: '15-Jun-10',
        37: '31-Jul-10',
        38: '15-Sep-10',
        39: '31-Oct-10',
        40: '16-Dec-10',
        41: '31-Jan-11',
        42: '18-Mar-11',
        43: '03-May-11'
    }

    def __init__(self, filename):
        """ Inspect a single file and get some basic info """
        super(SARAsset, self).__init__(filename)

        bname = os.path.basename(filename)
        self.tile = bname[10:17]
        self.sensor = bname[-9:-8] + bname[-15:-12]

        datafiles = self.datafiles()
        for f in datafiles:
            if f[-3:] == 'hdr':
                hdrfile = f
            if f[-4:] == 'date':
                datefile = f
                rootname = f[:-5]

        # unique to SARData (TODO - is this still used later?)
        self.hdrfile = hdrfile

        # Check if inspecting a file in the repository
        path = os.path.dirname(filename)
        if self.Repository.path() in path:
            date = datetime.datetime.strptime(os.path.basename(path), self.Repository._datedir).date()
            dates = date
            #VerboseOut('Date from repository = '+str(dates),4)
        else:
            # extract header and date image
            tfile = tarfile.open(filename)
            tfile.extract(hdrfile, path)
            hdrfile = os.path.join(path, hdrfile)
            os.chmod(hdrfile, 0664)
            meta = self._meta(hdrfile)
            tfile.extract(datefile, path)
            datefile = os.path.join(path, datefile)
            os.chmod(datefile, 0664)
            # Write ENVI header for date image
            List2File(meta['envihdr'], datefile + '.hdr')
            dateimg = gippy.GeoImage(datefile)
            dateimg.SetNoData(0)
            datevals = numpy.unique(dateimg.Read())
            dates = [self._launchdate[self.sensor[0]] + datetime.timedelta(days=int(d)) for d in datevals if d != 0]
            if not dates:
                RemoveFiles([hdrfile, datefile], ['.hdr', '.aux.xml'])
                raise Exception('%s: no valid dates' % bname)
            date = min(dates)
            dateimg = None
            RemoveFiles([hdrfile, datefile], ['.hdr', '.aux.xml'])
            #VerboseOut('Date from image: %s' % str(date),3)
            # If year provided check
            #if fname[7] == 'Y' and fname[8:10] != '00':
            #    ydate = datetime.datetime.strptime(fname[8:10], '%y')
            #    if date.year != ydate.year:
            #        raise Exception('%s: Date %s outside of expected year (%s)' % (fname, str(date),str(ydate)))
            # If widebeam check cycle dates
            if bname[7] == 'C':
                cdate = datetime.datetime.strptime(self._cycledates[int(bname[8:10])], '%d-%b-%y').date()
                if not (cdate <= date <= (cdate + datetime.timedelta(days=45))):
                    raise Exception('%s: Date %s outside of cycle range (%s)' % (bname, str(date), str(cdate)))
            #VerboseOut('%s: inspect %s' % (fname,datetime.datetime.now()-start), 4)
        self.date = dates
        self.rootname = rootname

    @classmethod
    def _meta(cls, hdrfile):
        """ Get some metadata from header file """
        hdr = File2List(hdrfile)
        meta = {}
        meta['proj'] = (
            'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984", SPHEROID["WGS_1984",6378137.0,298.257223563]],' +
            'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
        meta['size'] = [int(hdr[23]), int(hdr[24])]
        lat = [float(hdr[12]), float(hdr[14])]
        lat = [min(lat), max(lat)]
        meta['lat'] = lat
        lon = [float(hdr[13]), float(hdr[15])]
        lon = [min(lon), max(lon)]
        meta['lon'] = lon
        meta['res'] = [(lon[1] - lon[0]) / (meta['size'][0] - 1), (lat[1] - lat[0]) / (meta['size'][1] - 1)]
        meta['envihdr'] = [
            'ENVI', 'samples = %s' % meta['size'][0], 'lines = %s' % meta['size'][1],
            'bands = 1', 'header offset = 0', 'file type = ENVI Standard', 'data type = 12',
            'interleave = bsq', 'sensor type = Unknown', 'byte order = 0',
            'coordinate system string = ' + meta['proj'],
            'data ignore value = 0',
            'map info = {Geographic Lat/Lon, 1, 1, %s, %s, %s, %s}'
            % (lon[0], lat[1], meta['res'][0], meta['res'][1])]
        meta['CF'] = float(hdr[21])
        return meta

    def extract(self, filenames=[]):
        """ Extract filenames from asset and create ENVI header files """
        files = super(SARAsset, self).extract(filenames)
        for f in files:
            if f[-3:] == 'hdr':
                meta = self._meta(f)
        datafiles = {}
        for f in files:
            bname = os.path.basename(f)
            if f[-3:] != 'hdr':
                bandname = bname[len(self.rootname) + 1:]
                envihdr = copy.deepcopy(meta['envihdr'])
                if bandname in ['mask', 'linci']:
                    envihdr[6] = 'data type = 1'
                envihdr.append('band names={%s}' % bandname)
                List2File(envihdr, f + '.hdr')
            else:
                bandname = 'hdr'
            datafiles[bandname] = f
        return datafiles


class SARData(Data):
    """ Assets and products for a tile and date """
    name = 'SAR'
    version = '0.9.0'
    Asset = SARAsset

    _pattern = '*'
    _products = {
        'sign': {
            'description': 'Sigma nought (radar backscatter coefficient)',
        },
        'linci': {
            'description': 'Incident angles',
        },
        'date': {
            'description': 'Day of year array',
        },
    }

    def meta(self):
        """ Get metadata from headerfile """
        files = self.assets[''].datafiles()
        for f in files:
            if f[-3:] == 'hdr':
                hdr_bname = f
        files = self.assets[''].extract(filenames=[hdr_bname])
        return self.Asset._meta(files['hdr'])

    def find_files(self):
        """ Search path for valid files """
        filenames = super(SARData, self).find_files()
        filenames[:] = [f for f in filenames if os.path.splitext(f)[1] != '.hdr']
        return filenames

    def process(self, *args, **kwargs):
        """ Make sure all products have been pre-processed """
        products = super(SARData, self).process(*args, **kwargs)
        if len(products) == 0:
            return

        sensor = self.sensor_set[0]
        self.basename = self.basename + '_' + sensor
        # extract all data from archive
        datafiles = self.assets[''].extract()
        meta = self.meta()
        for key, val in products.requested.items():
            fname = os.path.join(self.path, self.basename + '_' + key)
            if val[0] == 'sign':
                bands = [datafiles[b] for b in ["sl_HH", "sl_HV"] if b in datafiles]
                img = gippy.GeoImage(bands)
                img.SetNoData(0)
                mask = gippy.GeoImage(datafiles['mask'], False)
                img.AddMask(mask[0] == 255)
                # apply date mask
                dateimg = gippy.GeoImage(datafiles['date'], False)
                dateday = (self.date - SARAsset._launchdate[sensor[0]]).days
                img.AddMask(dateimg[0] == dateday)
                #imgout = gippy.SigmaNought(img, fname, meta['CF'])
                imgout = gippy.GeoImage(fname, img, gippy.GDT_Float32)
                imgout.SetNoData(-32768)
                for b in range(0, imgout.NumBands()):
                    imgout.SetBandName(img[b].Description(), b + 1)
                    (img[b].pow(2).log10() * 10 + meta['CF']).Process(imgout[b])
                fname = imgout.Filename()
                img = None
                imgout = None
            if val[0] == 'linci':
                # Note the linci product DOES NOT mask by date
                os.rename(datafiles['linci'], fname)
                os.rename(datafiles['linci'] + '.hdr', fname + '.hdr')
            if val[0] == 'date':
                # Note the date product DOES NOT mask by date
                os.rename(datafiles['date'], fname)
                os.rename(datafiles['date'] + '.hdr', fname + '.hdr')
            self.AddFile(sensor, key, fname)
        # Remove unused files
        # TODO - checking key rather than val[0] (the full product suffix)
        if 'hdr' in datafiles:
            del datafiles['hdr']
        RemoveFiles(datafiles.values(), ['.hdr', '.aux.xml'])
