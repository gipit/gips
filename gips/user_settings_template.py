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

# GIPS user settings file

# this block will try to read in an environment (e.g., system-level or virtual)
# settings file so that individual settings can be changed
try:
    import gips.settings
    execfile(gips.settings.__file__.rstrip('c'))
except:
    raise Exception('There are no environment-level GIPS settings!')


# change email used for FTP
#EMAIL = ''


# to add in a database add a new key to the DATABASES dictionary
#DATABASES['mydatabase'] = {
#        'NAME': '',
#        'USER': '',
#        'PASSWORD': '',
#        'HOST': '',
#        'PORT': '5432',
#}


"""
# to add repository add new key to the REPOS dictionary
REPOS['dataname'] = {
   	# path to driver location (default to gips/data/dataname)
   	'driver': '',
   	# path to top level directory of data
  	'repopath': '',
   	# override location of tiles vector (defaults to gips/data/dataname/tiles.shp)
   'tiles': '',
   	#'tiles': 'mydatabase:mydatatype_tiles',		# database format
    #'tiles': '~/randomdir/dataname_tiles.shp'  	# file format

   # 'attribute name holding tileid in tiles vector'
   'tileid_attribute': '',
}
"""