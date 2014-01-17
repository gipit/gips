#!/usr/bin/env python

"""
setup for GIP and gippy
"""

#from distutils.core import setup, Extension
from setuptools import setup, Extension
import os, glob, filecmp, shutil
#os.environ['CC'] = 'g++';
#os.environ['CXX'] = 'g++';
#os.environ['CPP'] = 'g++';
#os.environ['LDSHARED'] = 'g++';

# Install GIP files
def install_gip(install):
    src = glob.glob('GIP/giputils/bin/Release/*')
    dst = [os.path.join('/usr/local/bin',os.path.basename(b)) for b in src]
    src.append('GIP/gip/bin/Release/libgip.so')
    dst.append('/usr/lib/libgip.so')
    for b in zip(src,dst):
        if not os.path.exists(b[1]) or not filecmp.cmp(b[0],b[1]): shutil.copy(b[0],b[1])

gippy_module = Extension(name = '_gippylib',
                    compiler='g++',
                    sources=['gippy/gippylib.i'],
                    swig_opts=['-c++', '-w509','-IGIP/gip'],
                    include_dirs=['GIP/gip'],
                    libraries=['gip','gdal','boost_system','boost_filesystem','X11'],
                    library_dirs=['GIP/gip/bin/Release'], #'/usr/lib','/usr/local/lib'],
                    extra_compile_args=['-fPIC'],
                    #extra_compile_args=['-fPIC -std=c++0x'],
                    ) 

setup (name = 'gippy',
    version = '1.0',
    description='Python bindings for GIP library',
    author='Matthew Hanson',
    author_email='mhanson@appliedgeosolutions.com',
    ext_modules = [gippy_module],
    py_modules = ['gippy.gippylib','gippy.atmosphere','gippy.GeoVector','gippy.gipit'],
    packages = ['gippy.algorithms','gippy.data'],
    dependency_links = ['https://github.com/matthewhanson/Py6S.git'],
    install_requires = ['Py6S','shapely==1.2.18'],
    entry_points = {
        'console_scripts': [
            'gipit = gippy.gipit:main',
            # Data scripts TODO - auto find and add
            'landsat = gippy.data.landsat:main',
            'cdl = gippy.data.cdl:main',
        ],
    },
    #cmdclass = {"install_gip": install_gip},
)
