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
setup for GIPS
"""

import os
from setuptools import setup
import gips.settings as settings
import glob
from gips.version import __version__

# console scripts
console_scripts = []
# Data scripts
for repo, cfg in settings.REPOS.items():
    if cfg['rootpath'] != '':
        console_scripts.append('%s = gips.data.%s:main' % (repo, repo.lower()))
# Algorithms
for f in glob.glob('gips/algorithms/*.py'):
    try:
        name = os.path.splitext(os.path.basename(f))[0]
        if name != '__init__':
            script = 'gips_%s = gips.algorithms.%s:main' % (name, name.lower())
            console_scripts.append(script)
    except:
        pass

scripts = []
if os.path.exists('bin'):
    files = glob.glob('bin/*.py')
    for f in files:
        scripts.append(f)

setup(
    name='gips',
    version=__version__,
    description='Geospatial Image Processing System',
    author='Matthew Hanson',
    author_email='mhanson@appliedgeosolutions.com',
    packages=['gips', 'gips.data', 'gips.algorithms'],
    install_requires=['Py6S>=1.5.0', 'shapely', 'gippy>=0.9.8', 'python-dateutil'],
    scripts=scripts,
    entry_points={'console_scripts': console_scripts},
)
