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

# Site files and data tiles vectors can be retrieved from a database
DATABASES = {
#    'tiles': {
#        'NAME': '',
#        'USER': '',
#        'PASSWORD': '',
#        'HOST': '',
#        'PORT': '5432',
#    }
}

# Used for anonymous FTP
EMAIL = '$EMAIL'


"""
# to add a repository add a new key to the REPOS dictionary, available params shown below

    'dataname': {
        # path to driver location (default to gips/data/dataname)
        'driver': '',
        # path to top level directory of data
        'repopath': '',
        # override location of tiles vector (defaults to gips/data/dataname/tiles.shp)
       'tiles': '',
        #'tiles': 'mydatabase:mydatatype_tiles',        # database format
        #'tiles': '~/randomdir/dataname_tiles.shp'      # file format

       # 'attribute name holding tileid in tiles vector'
       'tileid_attribute': '',
    }
"""


REPOS = {
    'AOD': {
        'rootpath': '$TLD/aod',
    },
    'Landsat': {
        'rootpath': '$TLD/landsat',
        'tile_attribute': 'pr',
        # Atmospheric correction
        '6S': False,
        'MODTRAN': False,
        # extract the files from tar.gz before processing instead of accessing directly
        'extract': False,
    },
    'Modis': {
        'rootpath': '$TLD/modis',
    },
    # these drivers tend to more specialized and experimental so turned off by default
    #'CDL': {
    #    'rootpath': '$TLD/cdl',
    #    'tile_attribute': '',
    #},
    #'SAR': {
    #    'rootpath': '$TLD/sar',
    #},
    #'SARAnnual': {
    #    'rootpath': '$TLD/sarannual',
    #},
    #'Merra': {
    #    'rootpath': '$TLD/Merra',
    #    'tile_attribute': 'tileid',
    #},
    #'Daymet': {
    #    'rootpath': '$TLD/daymet',
    #},
}
