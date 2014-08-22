#!/usr/bin/python

import argparse
import numpy as np

import gippy
from gips.algorithms.core import Algorithm
from gips.utils import VerboseOut


class Test(Algorithm):
    name = 'Test'
    __version__ = '0.1.0'

    def run(self, **kwargs):
        #arr = np.array([[1, 2], [3, 4], [5, 6]])
        arr = np.array([1, 2, 3, 4, 5, 6])

        VerboseOut('Original Array')
        print arr, arr.dtype, '\n'

        VerboseOut('Array w/o casting')
        outarr = gippy.test(arr)
        print outarr, outarr.dtype, '\n'

        dtypes = ['uint8', 'uint16', 'int16', 'uint32', 'int32', 'uint64', 'int64', 'float32', 'float64']
        for dt in dtypes:
            VerboseOut('Array as %s' % dt)
            outarr = gippy.test(arr.astype(dt))
            print outarr, outarr.dtype, '\n'

    @classmethod
    def parser(cls):
        parser = argparse.ArgumentParser(add_help=False)
        return parser


def main():
    Test.main()
