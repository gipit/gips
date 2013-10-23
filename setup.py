#!/usr/bin/env python

"""
setup for GIP and gippy
"""

from distutils.core import setup, Extension

gippy_module = Extension(name = '_gippylib',
                    sources=['gippy/gippylib.i'],
                    swig_opts=['-c++', '-w509','-IGIP/gip'],
                    include_dirs=['GIP/gip'],
                    libraries=['gip','gdal','boost_system','boost_program_options','boost_filesystem'],
                    extra_compile_args=['-fPIC'],
                    ) 

setup (name = 'gippy',
        version = '1.0',
        description='Python bindings for GIP library',
        author='Matthew Hanson',
        author_email='mhanson@appliedgeosolutions.com',
        ext_modules = [gippy_module],
        #py_modules = ['gippy.atmosphere','gippy.GeoVector'],
        #packages = ['gippy.algorithms','gippy.data']
        )
