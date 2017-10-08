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
from gips.utils import Colors, VerboseOut, import_data_class
from gips import utils
from gips.inventory import orm, dbinv


def main():
    title = Colors.BOLD + 'GIPS Data Archive Utility (v%s)' % gipsversion + Colors.OFF

    # argument parsing
    parser = GIPSParser(description=title)
    group = parser.add_argument_group('archive options')
    group.add_argument('--keep', help='Keep files after adding to archive', default=False, action='store_true')
    group.add_argument('--recursive', help='Iterate through subdirectories', default=False, action='store_true')
    group.add_argument(
        '--update',
        help='Update asset if newer version available, (must call gips_process to regenerate products',
        default=False,
        action='store_true'
    )
    args = parser.parse_args()

    utils.gips_script_setup(None, args.stop_on_error)

    with utils.error_handler('Data archive error'):
        print title
        cls = import_data_class(args.command)
        orm.setup() # set up DB orm in case it's needed for Asset.archive()
        # TODO archive accepts limited args, pass them in explicitly
        archived_assets = cls.Asset.archive(**vars(args))

        # if DB inventory is enabled, update it to contain the newly archived assets
        if orm.use_orm():
            for a in archived_assets:
                dbinv.update_or_add_asset(asset=a.asset, sensor=a.sensor, tile=a.tile, date=a.date,
                                          name=a.archived_filename, driver=cls.name.lower())

    utils.gips_exit()


if __name__ == "__main__":
    main()
