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

"""
setup for GIPS
"""

import os
from setuptools import setup
import shutil
import glob
import traceback
import imp

__version__ = imp.load_source('gips.version', 'gips/version.py').__version__

cfgpth = '/etc/gips'

# create cfgpath
try:
    if not os.path.exists(cfgpth):
        os.mkdir(cfgpth)
except OSError:
    # perhaps due to not root permissions but this may be a virtualenv so forge on ahead
    pass

# copies the GIPPY configuration file
try:
    configfile = os.path.join(cfgpth, 'settings.py')
    configtemplate = 'gips/settings.template.py'
    if not os.path.exists(configfile):
        shutil.copyfile(configtemplate, configfile)
except OSError:
    # perhaps due to not root permissions but this may be a virtualenv so forge on ahead
    pass

# copy the tile vectors
try:
    for d in glob.glob('data/*'):
        target = '/etc/gips/%s' % os.path.basename(d)
        if os.path.isdir(d) and not os.path.exists(target):
            shutil.copytree(d, target)
except:
    # perhaps due to not root permissions but this may be a virtualenv so forge on ahead
    pass


# collect console scripts
console_scripts = []
for f in glob.glob('gips/scripts/*.py'):
    try:
        name = os.path.splitext(os.path.basename(f))[0]
        if name not in ['__init__', 'core']:
            script = 'gips_%s = gips.scripts.%s:main' % (name, name.lower())
            console_scripts.append(script)
    except:
        print traceback.format_exc()

setup(
    name='gips',
    version=__version__,
    description='Geospatial Image Processing System',
    author='Matthew Hanson',
    author_email='matt.a.hanson@gmail.com',
    packages=['gips', 'gips.data', 'gips.scripts'],
    package_data={'': ['settings*py']},
    install_requires=['Py6S>=1.5.0', 'shapely', 'gippy', 'python-dateutil', 'pydap', 'scipy'],
    entry_points={'console_scripts': console_scripts},
)
