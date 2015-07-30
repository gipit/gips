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

# GIPS complete settings file

# Used for anonymous FTP
EMAIL = '$EMAIL'


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


REPOS = {
    'aod': {
        'repository': '$TLD/aod',
    },
    'landsat': {
        'repository': '$TLD/landsat',
        # Landsat specific settings
        '6S': False,            # atm correction for VIS/NIR/SWIR bands
        'MODTRAN': False,       # atm correction for LWIR
        'extract': False,       # extract files from tar.gz before processing instead of direct access
    },
    'modis': {
        'repository': '$TLD/modis',
    },
    # these drivers tend to more specialized and experimental so turned off by default
    #'cdl': {
    #    'repository': '$TLD/cdl',
    #},
    #'sar': {
    #    'repository': '$TLD/sar',
    #},
    #'sarannual': {
    #    'repository': '$TLD/sarannual',
    #},
    #'merra': {
    #    'repository': '$TLD/Merra',
    #},
    #'daymet': {
    #    'repository': '$TLD/daymet',
    #},
}


"""
# to add repository add new key to the REPOS dictionary
    'dataname': {
        # path to driver directory location (default to gips/data/dataname/ if not given)
        'driver': '',
        # path to top level directory of data
        'repository': '',
        # override location of tiles vector (default to gips/data/dataname/tiles.shp)
       'tiles': '',
        #'tiles': 'mydatabase:mydatatype_tiles',        # database format
        #'tiles': '~/randomdir/dataname_tiles.shp'      # file format
    }
"""