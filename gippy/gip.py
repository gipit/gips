#!/usr/bin/env python

import os, sys, argparse
import commands
import datetime

import shutil
import tempfile

import gippy
from gippy.GeoVector import GeoVector, NewShapefile, intersection

if __name__ == "__main__":
    prog = os.path.split(__file__)[1]
    dhf = argparse.ArgumentDefaultsHelpFormatter
    parser0 = argparse.ArgumentParser(prog=prog, description='Geospatial Image Processing', formatter_class=dhf)
    subparser = parser0.add_subparsers(dest='command')

    dbparser = argparse.ArgumentParser(add_help=False,formatter_class=dhf)
    group = dbparser.add_argument_group('Database connection options')
    group.add_argument('--dbname',help='Database name', default='geodata')
    group.add_argument('--dbhost',help='Database host', default='congo')
    group.add_argument('--dbuser',help='Database user', default='ags')
    group.add_argument('--dbpass',help='Database password', default='')
    group.add_argument('--dbport',help='Database port', default='5432')
    #group.add_argument('--dbtable',help='Database table name')

    parser = subparser.add_parser('inventory',parents=[dbparser],help='Show remote sensing data inventory',formatter_class=dhf)
    parser.add_argument('-s','--site',help='Vector file for region of interest', required=True)
    parser.add_argument('-d','--dates',help='Range of dates (YYYY-MM-DD,YYYY-MM-DD)')
    parser.add_argument('--days',help='Include only those that fall within these days of year (doy1,doy2)',default=None)

    parser = subparser.add_parser('create',help='Create site vector file',parents=[dbparser], formatter_class=dhf)
    parser.add_argument('layer', help='Data layer to use')
    parser.add_argument('-c','--col', help='Data layer to use', default='name')
    parser.add_argument('-v','--value',nargs='+',help='Criteria (col == value)', default='*')
    parser.add_argument('-w','--where', help='Where clause', default='')
    parser.add_argument('-o','--output',help='Output filename', default='site.shp')
    parser.add_argument('--overwrite',help='Overwrite output files', default=False, action='store_true')

    parser = subparser.add_parser('crop',help='Crop raster files to vector extents')
    parser.add_argument('files', nargs='*', help='Imagery files to extract site area from')
    parser.add_argument('-v','--vector',help='Shapefile of site to be extracted from imagery files', required=True)
    parser.add_argument('-r','--resolution', help='Resolution of output', type=float, default=None)
    #parser.add_argument('-w','--warp', help='Reproject raster files to match vector', action='store_true', default=False)
    parser.add_argument('-o', '--overwrite', help='Overwrite existing files', action='store_true', default=False)
    parser.add_argument('-a', '--alpha', help='Add alpha channel', action='store_true', default=False)

    parser = subparser.add_parser('int')
   
    args = parser0.parse_args()

    sensors = ['landsat',]
    for s in sensors: exec('import gippy.%s as %s' % (s,s))

    dbstr  = "PG:dbname=%s host=%s port=%s user=%s password=%s" % (args.dbname, args.dbhost, args.dbport, args.dbuser, args.dbpass)

    if args.command == 'inventory':
        for s in sensors:
            tiles = 
            tile_vector = GeoVector(dbstr, layer=eval('%s.tiles_vector' % s))
            #tiles = intersection(GeoVector(args.site),GeoVector(dbstr,layer=eval('%s.tile_vector' % s))
            #for tile,cover in tiles.items():
            #    if len(tile) == 5: tile = '0' + tile
            #    print "%s: %4.0f%%" % (key, val * 100)
            import pdb
            pdb.set_trace()
    elif args.command == 'crop':
        for f in args.files:
            start = datetime.datetime.now()
            fbase,ext = os.path.splitext(os.path.basename(f))
            fout = fbase + '_' + os.path.splitext(os.path.basename(args.vector))[0] + ext
            # update to get no data value from source....for that matter do this whole thing with API not cmdline utils
            cmd = 'gdalwarp %s %s -cutline %s -crop_to_cutline -srcnodata 0 -dstnodata -32768' % (f,fout,args.vector)
            if args.resolution: cmd = cmd+' -tr %s %s' % (args.resolution,args.resolution)
            if args.overwrite: cmd = cmd + ' -overwrite'
            if args.alpha: cmd = cmd + ' -dstalpha'
            out = commands.getstatusoutput(cmd)
            print '%s -> %s: %s' % (os.path.basename(f), os.path.basename(fout), datetime.datetime.now()-start)
    elif args.command == 'create':
        start = datetime.datetime.now()
        dbstr = "dbname=%s host=%s port=%s user=%s password=%s" % (args.dbname, args.dbhost, args.dbport, args.dbuser, args.dbpass)
        #sql = "select * FROM %s where %s.%s = '%s'" % (args.layer,args.layer,args.col,args.value[0])
        sql = "select * FROM %s where %s" % (args.layer,args.where)
        #if len(args.value) > 1: 
        #    for v in args.value[1:]: sql = sql + " OR %s.%s = '%s' " % (args.layer, args.col, v)
        #for a in args.admin1:
        #    sql = sql + ' admin1.name = %s ' % a
        #for a in args.admin2:
        #    sql = sql + ' usa_counties.name = %s ' % a
        cmd = 'ogr2ogr %s PG:"%s" -sql "%s"' % (args.output,dbstr,sql)
        if args.overwrite: cmd = cmd + ' -overwrite'
        print cmd
        out = commands.getstatusoutput(cmd)
        #print out
        print '%s: %s' % (args.output,datetime.datetime.now()-start)
