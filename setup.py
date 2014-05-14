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
import glob
from gipif.version import __version__

print 'GIPIF setup'

# console scripts
console_scripts = []
for repo, cfg in settings.REPOS.items():
    if cfg['rootpath'] != '':
        console_scripts.append('%s = gipif.data.%s:main' % (repo, repo.lower()))

scripts = []
if os.path.exists('bin'):
    files = glob.glob('bin/*.py')
    for f in files:
        scripts.append(f)

print 'scripts',scripts 

setup(
    name='gipif',
    version=__version__,
    description='Geospatial Image Processing and Inventory Framework',
    author='Matthew Hanson',
    author_email='mhanson@appliedgeosolutions.com',
    packages=['gipif', 'gipif.data'],
    install_requires=['Py6S', 'shapely', 'gippy>=0.9.4', 'python-dateutil'],
    scripts=scripts,
    entry_points={'console_scripts': console_scripts},
)
