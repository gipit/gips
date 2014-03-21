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
setup for GIP and gippy
"""

import os
from setuptools import setup, Extension
from setuptools.command.install import install
from setuptools.command.develop import develop
import numpy


class GIPinstall(install):
    def run(self):
        os.system('cd GIP; make; cd ..')
        install.run(self)


class GIPdevelop(develop):
    def initialize_options(self):
        #self.runtime_library_dirs = [os.path.abspath('GIP/bin/Release')]
        develop.initialize_options(self)
        #set_trace()
        self.ext_modules = [gippydev_module]

    def run(self):
        os.system('cd GIP; make; cd ..')
        develop.run(self)

#libgip = Extension(
#    name='libgip',
#    sources=['GIP/Atmosphere.cpp', 'GIP/GeoAlgorithms.cpp', 'GIP/GeoData.cpp',
#             'GIP/GeoImage.cpp', 'GIP/GeoRaster.cpp', 'GIP/GeoVector.cpp'],
#    include_dirs=['GIP'],
#    extra_compile_args=['-std=c++0x', '-Wall', '-fexceptions', '-fPIC', '-O2']
#)

gippy_module = Extension(
    name='_gippylib',
    sources=['gippy/gippylib.i'],
    swig_opts=['-c++', '-w509', '-IGIP'],
    include_dirs=['GIP', numpy.get_include()],
    libraries=['gip', 'gdal', 'boost_system', 'boost_filesystem'],  # ,'X11'],
    library_dirs=['GIP/bin/Release'],  # '/usr/lib','/usr/local/lib'],
    extra_compile_args=['-fPIC'],  # , '-std=c++0x'],
    #extra_compile_args=['-fPIC -std=c++0x'],
)
gippydev_module = gippy_module
gippydev_module.runtime_library_dirs = [os.path.abspath('GIP/bin/Release')]

setup(
    name='gippy',
    version='1.0',
    description='Python bindings for GIP library',
    author='Matthew Hanson',
    author_email='mhanson@appliedgeosolutions.com',
    ext_modules=[gippy_module],
    packages=['gippy', 'gippy.algorithms', 'gippy.data'],
    py_modules=['gippy.gippylib', 'gippy.atmosphere', 'gippy.GeoVector', 'gippy.gipit'],
    dependency_links=['https://github.com/matthewhanson/Py6S.git'],
    #install_requires = ['Py6S','shapely==1.2.18'],
    install_requires=['Py6S', 'shapely'],
    data_files=[('/usr/lib', ['GIP/bin/Release/libgip.so'])],
    entry_points={
        'console_scripts': [
            'gipit = gippy.gipit:main',
            # Data scripts TODO - auto find and add
            'landsat = gippy.data.landsat:main',
            'CDL = gippy.data.cdl:main',
            'SAR = gippy.data.sar:main',
            'SARannual= gippy.data.sarannual:main',
            'modis = gippy.data.modis:main',
        ],
    },
    cmdclass={
        "develop": GIPdevelop,
        "install": GIPinstall
    }
)
