#!/usr/bin/env python
################################################################################
#    GIPS: Geospatial Image Processing System
#
#    AUTHOR: Matthew Hanson
#    EMAIL:  mhanson@ags.io
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
requirements = []
# Data scripts
for repo, cfg in settings.REPOS.items():
    if cfg['rootpath'] != '':
        console_scripts.append('%s = gips.data.%s:main' % (repo, repo.lower()))
        try:
            exec('from gips.data.%s import requirements as reqs' % repo.lower())
            requirements.extend(reqs)
        except:
            pass
# Algorithms
for f in glob.glob('gips/algorithms/*.py'):
    try:
        name = os.path.splitext(os.path.basename(f))[0]
        if name not in ['__init__', 'core']:
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
    author_email='mhanson@ags.io',
    packages=['gips', 'gips.data', 'gips.algorithms'],
    # need pydap for merra
    install_requires=requirements.extend(['shapely', 'gippy>=1.0.0', 'python-dateutil']),
    scripts=scripts,
    entry_points={'console_scripts': console_scripts},
)
