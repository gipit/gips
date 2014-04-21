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
from setuptools import setup
import gipif.settings as settings

from pdb import set_trace

# console scripts
console_scripts = []
for repo, cfg in settings.REPOS.items():
    if cfg['rootpath'] != '':
        console_scripts.append('%s = gipif.data.%s:main' % (repo, repo.lower()))

setup(
    name='gipif',
    version='1.0',
    description='Geospatial Image Processing and Inventory Framework',
    author='Matthew Hanson',
    author_email='mhanson@appliedgeosolutions.com',
    packages=['gipif'],
    #py_modules=['gippy.GeoVector', 'gippy.gipit'],
    #dependency_links=[
    #    'https://github.com/robintw/Py6S.git'
    #    'https://bitbucket.org/appliedgeosolutions/gippy.git'],
    #install_requires = ['Py6S','shapely==1.2.18'],
    install_requires=['Py6S', 'shapely', 'gippy'],
    entry_points={'console_scripts': console_scripts},
)
