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

#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os, glob, filecmp, shutil
import numpy
#os.environ['CC'] = 'g++';
#os.environ['CXX'] = 'g++';
#os.environ['CPP'] = 'g++';
#os.environ['LDSHARED'] = 'g++';

# Install GIP files
def install_gip(install):
    src = glob.glob('GIP/giputils/bin/Release/*')
    dst = [os.path.join('/usr/local/bin',os.path.basename(b)) for b in src]
    src.append('GIP/bin/Release/libgip.so')
    dst.append('/usr/lib/libgip.so')
    for b in zip(src,dst):
        if not os.path.exists(b[1]) or not filecmp.cmp(b[0],b[1]): shutil.copy(b[0],b[1])

libgip = Extension(name = 'libgip',
            sources=['GIP/Atmosphere.cpp','GIP/GeoAlgorithms.cpp','GIP/GeoData.cpp',
                'GIP/GeoImage.cpp','GIP/GeoRaster.cpp','GIP/GeoVector.cpp'],
            include_dirs=['GIP'],
            extra_compile_args=['-std=c++0x','-Wall','-fexceptions','-fPIC','-O2']
        )

gippy_module = Extension(name = '_gippylib',
                    #compiler='g++',
                    sources=['gippy/gippylib.i'],
                    swig_opts=['-c++', '-w509','-IGIP'],
                    include_dirs=['GIP',numpy.get_include()],
                    libraries=['gip','gdal','boost_system','boost_filesystem'], #,'X11'],
                    library_dirs=['GIP/bin/Release'], #'/usr/lib','/usr/local/lib'],
                    extra_compile_args=['-fPIC'], #, '-std=c++0x'],
                    #extra_compile_args=['-fPIC -std=c++0x'],
                    ) 

setup (name = 'gippy',
    version = '1.0',
    description='Python bindings for GIP library',
    author='Matthew Hanson',
    author_email='mhanson@appliedgeosolutions.com',
    ext_modules = [gippy_module],
    packages = ['gippy','gippy.algorithms','gippy.data'],
    py_modules = ['gippy.gippylib','gippy.atmosphere','gippy.GeoVector','gippy.gipit'],
    dependency_links = ['https://github.com/matthewhanson/Py6S.git'],
    #install_requires = ['Py6S','shapely==1.2.18'],
    install_requires = ['Py6S','shapely'],
    #data_files = [('/usr/local/lib',['GIP/bin/Release/libgip.so'])],
    entry_points = {
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
    #cmdclass = {"install_gip": install_gip},
)
