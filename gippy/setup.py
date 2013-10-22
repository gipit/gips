#!/usr/bin/env python

"""
setup for GIP and gippy
"""

from distutils.core import setup, Extension

gippy_module = Extension(name = '_gippy',
                    sources=['gippylib.i'],
                    swig_opts=['-c++', '-w509','-I../GIP/gip'],
                    include_dirs=['../GIP/gip'],
                    ) 

setup (name = 'gippy',
        version = '0.1',
        ext_modules = [gippy_module],
        py_modules = ["gippy"],
        )
