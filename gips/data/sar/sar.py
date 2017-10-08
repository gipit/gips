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
import re
import datetime
import tarfile
import copy
import numpy

import gippy
from gips.data.core import Repository, Asset, Data
from gips.utils import File2List, List2File, RemoveFiles
from gips import utils


class sarRepository(Repository):
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


class sarAsset(Asset):
    """ Single original file """
    Repository = sarRepository

    _launchdate = {
        'alos1': datetime.date(2006, 1, 24),
        'jers': datetime.date(1992, 2, 11),
        'alos2': datetime.date(2014, 5, 24),
    }

    _sensors = {
        'AFBS': {
            'description': 'PALSAR FineBeam Single Polarization',
            'startdate': _launchdate['alos1'],
        },
        'AFBD': {
            'description': 'PALSAR FineBeam Dual Polarization',
            'startdate': _launchdate['alos1'],
        },
        'AWB1': {
            'description': 'PALSAR WideBeam (ScanSAR Short Mode)',
            'startdate': _launchdate['alos1'],
        },
        'AWBD': {
            'description': ('PALSAR-2 WideBeam '
                            '(ScanSAR w/Ortho Slope Correction)'),
            'startdate': _launchdate['alos2'],
        },
        'JFBS': {
            'description': 'JERS-1 FineBeam Single Polarization',
            'startdate': _launchdate['jers'],
        },
    }

    # N.B.> JAXA spec incl. blank lines (i.e. 16)
    __JAXA_spec_hdr_lines = {
        'prod': 0,
        'user': 1,
        'obsv': 2,
        'tile': 3,
        'op_mode': 4,
        'path_type': 5,
        'satellite': 6,
        'pix_spacing_meters': 7,
        'direction': 8,
        'resamp_meth': 9,
        'slope_correct': 10,
        'ul_lat': 12,
        'ul_lon': 13,
        'lr_lat': 14,
        'lr_lon': 15,
        'proj': 17,
        'y_res_arcsec': 18,
        'x_res_arcsec': 19,
        'cal_factor': 21,
        'rows': 23,
        'cols': 24,
        'processing_date': 26,
        'processing_time': 27,
    }
    _assets = {
        'alos1': {
            # KC_017-C27N00E100WB1ORSA1.tar.gz     # old
            'pattern': (
                r'^KC_(?P<userid>[0-9]{3})-'
                r'(?P<year_or_cycle>[CY])'
                r'(?P<cyid>[0-9]+)'
                r'(?P<tile>[NS][0-9]{2}[EW][0-9]{3})'
                r'(?P<mode>[FWP][LB][1DSR])'
                r'(?P<pathtype>OR[SM])'
                r'(?P<satellite>[AJ])'
                r'(?P<serialno>[1-9])'
                r'\.tar\.gz$'
            ),
            'cycledates': {
                7: '20-Oct-06', 8: '05-Dec-06', 9: '20-Jan-07',
                10: '07-Mar-07', 11: '22-Apr-07', 12: '07-Jun-07',
                13: '23-Jul-07', 14: '07-Sep-07', 15: '23-Oct-07',
                16: '08-Dec-07', 17: '23-Jan-08', 18: '09-Mar-08',
                19: '24-Apr-08', 20: '09-Jun-08', 21: '25-Jul-08',
                22: '09-Sep-08', 23: '25-Oct-08', 24: '10-Dec-08',
                25: '25-Jan-09', 26: '12-Mar-09', 27: '27-Apr-09',
                28: '12-Jun-09', 29: '28-Jul-09', 30: '12-Sep-09',
                31: '28-Oct-09', 32: '13-Dec-09', 33: '28-Jan-10',
                34: '15-Mar-10', 35: '30-Apr-10', 36: '15-Jun-10',
                37: '31-Jul-10', 38: '15-Sep-10', 39: '31-Oct-10',
                40: '16-Dec-10', 41: '31-Jan-11', 42: '18-Mar-11',
                43: '03-May-11',
            },
            'hdr_lines': copy.deepcopy(__JAXA_spec_hdr_lines),
        },
        'alos2': {
            # KC_999-C045DRN00E115WBDORSA1.tar.gz  # new
            'pattern': (
                r'^KC_(?P<userid>[0-9]{3})-'
                r'(?P<year_or_cycle>C)'
                r'(?P<cyid>[0-9]+)'
                r'(?P<format>[AD])'
                r'(?P<lookdir>[LR])'
                r'(?P<tile>(?P<lat>[NS][0-9]{2})(?P<lon>[EW][0-9]{3}))'
                r'(?P<mode>WBD)'
                r'ORS'
                r'(?P<satellite>[AJ])'
                r'(?P<serialno>[0-9])'
                r'\.tar\.gz$'
            ),
            'hdr_lines': copy.deepcopy(__JAXA_spec_hdr_lines),
        },
    }
    _assets['alos2']['hdr_lines'].update(
        {
            # HACK due to JAXA assets ignoring spec and omitting blank lines
            # (20 is supposed to be blank, but ALOS2 assets seem to ignore this)
            'cal_factor': 20,
            'rows': 22,
            'cols': 23,
            'processing_date': 25,
            'processing_time': 26,
        }
    )
    _defaultresolution = [0.000834028356964, 0.000834028356964]

    def __init__(self, filename):
        """Inspect a single file and get some basic info.

        **Nota Bene**: This driver is not gips.datahandler compliant due to
                       requiring the actual file on hand in order to
                       instantiate a class. As a possible solution, we could
                       add a date argument to the constructor, since any time
                       one is constructing an asset w/o the file, they likely
                       know the date for which they remote file is specified.
        """
        super(sarAsset, self).__init__(filename)

        bname = os.path.basename(filename)
        mats = {}
        for a in self._assets:
            m = re.match(self._assets[a]['pattern'], bname)
            if m:
                mats[a] = m

        if not mats:
            raise Exception(
                "{} doesn't match asset naming convention".format(bname)
            )
        elif len(mats) > 1:
            raise Exception('{} matches pattern for: ' + ','.join(mats))

        self.asset, m = mats.items()[0]
        self.tile = m.group('tile')
        self.sensor = m.group('satellite') + m.group('mode')

        self._version = int(m.group('serialno'))

        self.is_cycle = m.group('year_or_cycle') == 'C'
        self.cyid = int(m.group('cyid'))
        self.bname = bname

        # Check if inspecting a file in the repository
        path = os.path.dirname(filename)
        self._meta_dict = None
        if self.Repository.path() in path:
            date = datetime.datetime.strptime(
                os.path.basename(path),
                self.Repository._datedir
            ).date()
            self.date = date
            #VerboseOut('Date from repository = '+str(dates),4)
        elif os.path.exists(filename):
            meta = self.get_meta_dict()
            self.date = meta['min_date']


    def get_meta_dict(self):
        if not self._meta_dict:
            self._proc_meta()
        assert self._meta_dict
        return copy.deepcopy(self._meta_dict)

    def _jaxa_opener(self, f, path=None):
        if type(f) in (str, unicode):
            f = (f,)
        paths = self.extract(f, path=path).values()
        img = gippy.GeoImage(paths)
        return img

    def _proc_meta(self):
        """ Get some metadata from header file """
        # add asset specific hdr file line number keys to local namespace
        l = self._assets[self.asset]['hdr_lines']
        ###############
        ###############
        datafiles = self.datafiles()
        for f in datafiles:
            if f[-3:] == 'hdr':
                hdrfile = f
            if f[-4:] == 'date':
                datefile = f
                self.rootname = f[:-5]
        hdrfile = super(sarAsset, self).extract((hdrfile,))[0]
        hdr = File2List(hdrfile)
        RemoveFiles((hdrfile,))
        ##### LOAD META
        meta = {}
        meta['proj'] = (
            'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984", SPHEROID["WGS_1984",6378137.0,298.257223563]],' +
            'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')

        meta['size'] = [int(hdr[l['rows']]), int(hdr[l['cols']])]
        lat = sorted([float(hdr[l['ul_lat']]), float(hdr[l['lr_lat']])])
        meta['lat'] = lat
        lon = sorted([float(hdr[l['ul_lon']]), float(hdr[l['lr_lon']])])
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
        meta['CF'] = float(hdr[l['cal_factor']])

        # N.B.: self._meta_dict isn't complete, but is complete enough to open
        #       via _jaxa_opener.  Only need to get date to complete the it.
        self._meta_dict = meta

        dateimg = self._jaxa_opener(datefile)
        dateimg.SetNoData(0)
        datevals = numpy.unique(dateimg.Read())
        dateimg = None
        RemoveFiles((datefile,), ['.hdr', '.aux.xml'])
        dates = [
            (self._sensors[self.sensor]['startdate'] +
             datetime.timedelta(days=int(d)))
            for d in datevals if d != 0
        ]
        if not dates:
            raise Exception('%s: no valid dates' % self.bname)
        date = min(dates)
        self._meta_dict['min_date'] = date
        #VerboseOut('Date from image: %s' % str(date),3)
        # If year provided check
        #if fname[7] == 'Y' and fname[8:10] != '00':
        #    ydate = datetime.datetime.strptime(fname[8:10], '%y')
        #    if date.year != ydate.year:
        #        raise Exception('%s: Date %s outside of expected year (%s)' % (fname, str(date),str(ydate)))
        # If widebeam check cycle dates
        if self.asset == 'alos1' and self.is_cycle:
            cycledates = self._assets[self.asset]['cycledates']
            cdate = datetime.datetime.strptime(
                cycledates[self.cyid], '%d-%b-%y'
            ).date()
            delta = (date - cdate).days
            utils.verbose_out(
                '{}: {} days different between datearray and cycledate'
                .format(self.bname, delta)
            )
            if not (0 <= delta <= 45):
                raise Exception('%s: Date %s outside of cycle range (%s)' % (self.bname, str(date), str(cdate)))
        #VerboseOut('%s: inspect %s' % (fname,datetime.datetime.now()-start), 4)


    def extract(self, filenames=[], path=None):
        """ Extract filenames from asset and create ENVI header files """
        files = super(sarAsset, self).extract(filenames, path=path)
        meta = self.get_meta_dict()
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


class sarData(Data):
    """ Assets and products for a tile and date """
    name = 'SAR'
    version = '0.9.0'
    Asset = sarAsset

    _products = {
        'sign': {
            'description': 'Sigma nought (radar backscatter coefficient)',
            'assets': ['alos1', 'alos2'],
        },
        'mask': {
            'description': ('Mask (0:nodata, 50:water, '
                            '100:layover, 150:shadow, 255:land)'),
            'assets': ['alos1', 'alos2'],
        },
        'linci': {
            'description': 'Incident angles',
            'assets': ['alos1', 'alos2'],
        },
        'date': {
            'description': 'Day of year array',
            'assets': ['alos1', 'alos2'],
        },
    }

    def find_files(self):
        """ Search path for valid files """
        filenames = super(sarData, self).find_files()
        filenames[:] = [f for f in filenames if os.path.splitext(f)[1] != '.hdr']
        return filenames

    @Data.proc_temp_dir_manager
    def process(self, *args, **kwargs):
        """
        Make sure all products have been pre-processed.
        
        TODO: add metadata to output images.
        """
        products = super(sarData, self).process(*args, **kwargs)
        if len(products) == 0:
            return

        sensor = self.sensor_set[0]    # for a given time-space, there should
        #                              # only be a look from one sensor,
        asset = self.assets.keys()[0]  # and one asset.

        datafiles = self.assets[asset].extract(path=self._temp_proc_dir)

        #datafiles = self.assets[asset].datafiles()
        self.basename = self.basename + '_' + sensor
        for key, val in products.requested.items():
            #fname = os.path.join(self.path, self.basename + '_' + key)
            fname = self.temp_product_filename(sensor, key)
            jo = lambda fps: self.assets[asset]._jaxa_opener(
                    fps, path=self._temp_proc_dir)

            if val[0] == 'mask':
                mask = jo(datafiles['mask'])
                imgout = mask.Process(fname)
                imgout[0].SetNoData(0)
                imgout.Process()
            elif val[0] == 'sign':
                # extract all data from archive
                bands = [datafiles[b] for b in ["sl_HH", "sl_HV"] if b in datafiles]
                # bands = [b for b in datafiles
                #          if any((b.endswith(sl) for sl in ["sl_HH", "sl_HV"]))
                # ]
                img = jo(bands)
                img.SetNoData(0)
                mask = jo(datafiles['mask'])
                mask[0] = mask[0].BXOR(150.) > 0
                img.AddMask(mask[0])
                # apply date mask
                dateimg = jo(datafiles['date'])
                dateday = (
                    self.date -
                    self.assets[asset]._sensors[sensor]['startdate']).days
                img.AddMask(dateimg[0] == dateday)
                imgout = gippy.GeoImage(fname, img, gippy.GDT_Float32)
                imgout.SetNoData(-32768)
                for b in range(0, imgout.NumBands()):
                    imgout.SetBandName(img[b].Description(), b + 1)
                    (
                        img[b].pow(2).log10() * 10 +
                        self.assets[asset].get_meta_dict()['CF']
                    ).Process(imgout[b])
                fname = imgout.Filename()
            elif val[0] == 'linci':
                # Note the linci product DOES NOT mask by date
                img = gippy.GeoImage(datafiles['linci'])
                img.Process(fname)
            elif val[0] == 'date':
                # Note the date product DOES NOT mask by date
                img = gippy.GeoImage(datafiles['date'])
                img.Process(fname)
            else:
                raise Exception('Unrecognized product: ' + key)
            archive_fp = self.archive_temp_path(fname)
            self.AddFile(sensor, key, archive_fp)
        # Remove unused files
        # TODO - checking key rather than val[0] (the full product suffix)
        if 'hdr' in datafiles:
            del datafiles['hdr']
