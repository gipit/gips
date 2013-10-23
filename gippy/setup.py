#!/usr/bin/env python

"""
setup for GIP and gippy
"""

from distutils.core import setup, Extension



gippy_module = Extension(name = 'gippy',
                    sources=['gippylib.i'],
                    swig_opts=['-c++', '-w509','-I../GIP/gip'],
                    include_dirs=['../GIP/gip'],
                    ) 

setup (name = 'gippy',
        version = '1.0',
        description='Python bindings for GIP library',
        author='Matthew Hanson',
        author_email='mhanson@appliedgeosolutions.com',
        ext_modules = [gippy_module],
        py_modules = ['gippy.atmosphere','gippy.GeoVector'],
        packages = ['gippy.algorithms','gippy.data']
        )
