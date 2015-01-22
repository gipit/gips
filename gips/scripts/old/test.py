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

import gips
from gips.algorithms.core import Algorithm
from gips.utils import VerboseOut
import traceback


class Test(Algorithm):
    name = 'Test'
    __version__ = '0.1.0'

    def run(self, repos=None, **kwargs):
        repos = repos if repos is not None else gips.settings.REPOS
        for repo in repos:
            try:
                exec('from gips.data.%s import test' % repo.lower())
                print
                exec('test()')
            except Exception, e:
                VerboseOut('\n%s error: %s' % (repo, e))
                VerboseOut(traceback.format_exc(), 3)
                #print '\n%s: no test shapefile' % repo
                pass

    @classmethod
    def parser(cls, parser0):
        parser0.add_argument('--repos', help='Repositories to test', nargs='*', default=None)
        return parser0


def main():
    Test.main()
