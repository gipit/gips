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

"""
setup for GIPIF
"""

import os
import shutil
from setuptools import setup, Extension
from setuptools.command.install import install
from setuptools.command.develop import develop
from copy import deepcopy
import numpy


setup(
    name='gipif',
    version='1.0',
    description='Geospatial Image Processing and Inventory Framework',
    author='Matthew Hanson',
    author_email='mhanson@appliedgeosolutions.com',
    packages=['gipif'],
    #py_modules=['gippy.GeoVector', 'gippy.gipit'],
    dependency_links=[
        'https://github.com/robintw/Py6S.git'
        'https://bitbucket.org/appliedgeosolutions/gippy.git'],
    #install_requires = ['Py6S','shapely==1.2.18'],
    install_requires=['Py6S', 'shapely', 'gippy'],
    entry_points={
        'console_scripts': [
            'landsat = gipif.data.landsat:main',
            'CDL = gipif.data.cdl:main',
            'SAR = gip.data.sar:main',
            'SARannual= gipif.data.sarannual:main',
            'modis = gipif.data.modis:main',
            'atmos = gipif.data.atmos:main',
        ],
    },
)
