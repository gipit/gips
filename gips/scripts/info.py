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

from gips import __version__ as gipsversion
from gips.parsers import GIPSParser
from gips.data.core import data_class
from gips.utils import Colors, VerboseOut


def main():
    title = Colors.BOLD + 'GIPS Data Repositories (v%s)' % (gipsversion) + Colors.OFF

    # argument parsing
    parser = GIPSParser(description=title)
    args = parser.parse_args()

    try:
        print title
        cls = data_class(args.command)
        cls.print_products()
    except Exception, e:
        import traceback
        VerboseOut(traceback.format_exc(), 4)
        print 'GIPS Data Repository error: %s' % e


if __name__ == "__main__":
    main()
