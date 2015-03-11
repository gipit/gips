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

DATABASES = {
    'tiles': {
        'NAME': '',
        'USER': '',
        'PASSWORD': '',
        'HOST': '',
        'PORT': '5432',
    }
}

# Used for anonymous FTP
EMAIL = ''

REPOS = {
    'AOD': {
        'rootpath': '/data/aod',
    },
    'CDL': {
        'rootpath': '',
	    'tile_attribute': '',
    },
    'Landsat': {
        'rootpath': '/data/landsat',
	    'tile_attribute': 'pr',
        # Atmospheric correction
        '6S': False,
        'MODTRAN': False,
        # Extract the files from tar.gz before processing (opposed to accessing directly)
        'extract': False,
    },
    'Modis': {
        'rootpath': '/data/modis',
    },
    'SAR': {
        'rootpath': '',
    },
    'SARAnnual': {
        'rootpath': '',
    },
    'Merra': {
        'rootpath': '',
	    'tile_attribute': 'tileid',
    },
    'Daymet': {
        'rootpath': '',
    },
}
