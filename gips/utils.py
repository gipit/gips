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
from __future__ import print_function

import imp
import sys
import os
import re
import errno
from contextlib import contextmanager
import tempfile
import commands
import shutil
import traceback
import datetime

import gippy
from gippy import GeoVector


class Colors():
    _c = '\033['
    OFF     = _c + '0m'
    # Font styles
    BOLD    = _c + '1m'
    UNDER   = _c + '4m'
    REV     = _c + '7m'
    # Text colors
    BLACK   = _c + '30m'
    RED     = _c + '31m'
    GREEN   = _c + '32m'
    YELLOW  = _c + '33m'
    BLUE    = _c + '34m'
    PURPLE  = _c + '35m'
    CYAN    = _c + '36m'
    WHITE   = _c + '37m'
    # Background colors
    _BLACK  = _c + '40m'
    _RED    = _c + '41m'
    _GREEN  = _c + '42m'
    _YELLOW = _c + '43m'
    _BLUE   = _c + '44m'
    _PURPLE = _c + '45m'
    _CYAN   = _c + '46m'
    _WHITE  = _c + '47m'


def verbose_out(obj, level=1, stream=sys.stdout):
    """print(obj) but only if the user's chosen verbosity level warrants it.

    Print to stdout by default, but select any stream the user wishes.  Finally
    if the obj is a list or tuple, print each contained object consecutively on
    separate lines.
    """
    #TODO: Add real documentation of rules regarding levels used within
    #      GIPS. Levels 1-4 are used frequently.  Setting `-v5` is
    #      "let me see everything" level.
    if gippy.Options.Verbose() >= level:
        if not isinstance(obj, (list, tuple)):
            obj = [obj]
        for o in obj:
            print(o, file=stream)

VerboseOut = verbose_out # VerboseOut name is deprecated

##############################################################################
# Filesystem functions
##############################################################################

def File2List(filename):
    """Return contents of file as a list of lines, sans newlines."""
    f = open(filename)
    txt = f.readlines()
    txt2 = []
    for t in txt:
        txt2.append(t.rstrip('\n'))
    return txt2


def List2File(lst, filename):
    """Overwrite the given file with the contents of the list.

    Each item in the list is given a trailing newline.
    """
    f = open(filename, 'w')
    f.write('\n'.join(lst) + '\n')
    f.close()


def remove_files(filenames, extensions=()):
    """Remove the given files and all permutations with the given extensions.

    So remove_files(['a.hdf', 'b.hdf'], ['.index', '.aux.xml']) attempts to
    these files:  a.hdf, b.hdf, a.hdf.index, a.hdf.aux.xml, b.hdf,
    b.hdf.index, and b.hdf.aux.xml.  Doesn't raise an error if any file doesn't exist.
    """
    for f in (list(filenames) + [f + e for f in filenames for e in extensions]):
        with error_handler(continuable=True, msg_prefix="Error removing '{}'".format(f)):
            if os.path.isfile(f):
                os.remove(f)

RemoveFiles = remove_files # RemoveFiles name is deprecated

def basename(str):
    """Return the input string, stripped of directories and extensions.

    So, basename('/home/al-haytham/book-of-optics.pdf') returns
    'book-of-optics'.
    """
    return os.path.splitext(os.path.basename(str))[0]


def mkdir(dname):
    """ Create a directory if it doesn't exist """
    if not os.path.exists(dname):
        os.makedirs(dname)
    return dname


def link(src, dst, hard=False):
    """ Create link in this directory """
    if os.path.lexists(dst):
        os.remove(dst)
    # link path path relative to dst
    if hard:
        os.link(src, os.path.abspath(dst))
    else:
        os.symlink(os.path.relpath(src, os.path.dirname(dst)), os.path.abspath(dst))
    return dst

@contextmanager
def make_temp_dir(suffix='', prefix='tmp', dir=None):
    """Context manager to create then delete a temporary directory.

    Arguments are the same as tempfile.mkdtemp, which it calls.  Yields
    the absolute pathname to the new directory.  Deletes the directory
    at the exit of the context, regardless of exceptions raised in the
    context.
    """
    absolute_pathname = tempfile.mkdtemp(suffix, prefix, dir)
    try:
        yield absolute_pathname
    finally:
        shutil.rmtree(absolute_pathname)


def find_files(regex, path='.'):
    """Find filenames in the given directory that match the regex.

    Returns a list of matching filenames; each includes the given path.
    Only regular files and symbolic links to regular files are returned.
    """
    compiled_re = re.compile(regex)
    ret = []
    for f in os.listdir(path):
        fpath = os.path.join(path, f)
        if os.path.isfile(fpath) and compiled_re.match(f):
            ret.append(fpath)
    return ret


##############################################################################
# Settings functions
##############################################################################


def settings():
    """ Retrieve GIPS settings - first from user, then from system """
    settings_path = os.path.expanduser('~/.gips/settings.py')
    if os.path.isfile(settings_path):
        with error_handler("Error loading '{}'".format(settings_path)):
            # import user settings first
            src = imp.load_source('settings', settings_path)
            return src
    with error_handler("gips.settings not found; consider running gips_config"):
        import gips.settings
        return gips.settings


def create_environment_settings(repos_path, email=''):
    """ Create settings file and data directory """
    from gips.settings_template import __file__ as src
    cfgpath = os.path.dirname(__file__)
    cfgfile = os.path.join(cfgpath, 'settings.py')
    if src[-1] == 'c':
        src = src[:-1]
    if os.path.exists(cfgfile):
        return cfgfile
    with open(cfgfile, 'w') as fout:
        with open(src, 'r') as fin:
            for line in fin:
                fout.write(line.replace('$TLD', repos_path).replace('$EMAIL', email))
    return cfgfile


def create_user_settings(email=''):
    """ Create a settings file using the included template and the provided top level directory """
    from gips.user_settings_template import __file__ as src
    if src[-1] == 'c':
        src = src[:-1]
    cfgpath = os.path.expanduser('~/.gips')
    if not os.path.exists(cfgpath):
        os.mkdir(cfgpath)
    cfgfile = os.path.join(cfgpath, 'settings.py')
    if os.path.exists(cfgfile):
        raise Exception('User settings file already exists: %s' % cfgfile)
    with open(cfgfile, 'w') as fout:
        with open(src, 'r') as fin:
            for line in fin:
                fout.write(line)
    return cfgfile


def create_repos():
    """ Create any necessary repository directories """
    repos = settings().REPOS
    for key in repos.keys():
        repo = import_repository_class(key)
        for d in repo._subdirs:
            path = os.path.join(repos[key]['repository'], d)
            if not os.path.isdir(path):
                os.makedirs(path)


def data_sources():
    """ Get enabled data sources (and verify) from settings """
    sources = {}
    repos = settings().REPOS
    for key in sorted(repos.keys()):
        if not os.path.isdir(repos[key]['repository']):
            raise Exception('ERROR: archive %s is not a directory or is not available' % key)
        with error_handler(continuable=True):
            repo = import_repository_class(key)
            sources[key] = repo.description
    return sources


def import_data_module(clsname):
    """ Import a data driver by name and return as module """
    import imp
    path = settings().REPOS[clsname].get('driver', '')
    if path == '':
        path = os.path.join( os.path.dirname(__file__), 'data', clsname)
    with error_handler('Error loading driver ' + clsname):
        fmtup = imp.find_module(clsname, [path])
        mod = imp.load_module(clsname, *fmtup)
        return mod


def import_repository_class(clsname):
    """ Get clsnameRepository class object """
    mod = import_data_module(clsname)
    exec('repo = mod.%sRepository' % clsname)
    return repo


def import_data_class(clsname):
    """ Get clsnameData class object """
    mod = import_data_module(clsname)
    exec('repo = mod.%sData' % clsname)
    # prevent use of database inventory for certain incompatible drivers
    from gips.inventory import orm
    orm.driver_for_dbinv_feature_toggle = repo.name.lower()
    return repo


##############################################################################
# Geospatial functions
##############################################################################

def open_vector(fname, key="", where=''):
    """Open vector or feature, returned as a gippy GeoVector or GeoFeature."""
    parts = fname.split(':')
    if len(parts) == 1:
        vector = GeoVector(fname)
        vector.SetPrimaryKey(key)
    else:
        # or it is a database
        if parts[0] not in settings().DATABASES.keys():
            raise Exception("%s is not a valid database" % parts[0])
        db = settings().DATABASES[parts[0]]
        filename = ("PG:dbname=%s host=%s port=%s user=%s password=%s" %
                    (db['NAME'], db['HOST'], db['PORT'], db['USER'], db['PASSWORD']))
        vector = GeoVector(filename, parts[1])
        vector.SetPrimaryKey(key)
    if where != '':
        # return array of features
        return vector.where(where)
    else:
        return vector

from shapely.wkt import loads as wktloads
from osr import SpatialReference, CoordinateTransformation
from ogr import CreateGeometryFromWkt


def transform_shape(shape, ssrs, tsrs):
    """ Transform shape from ssrs to tsrs (all wkt) and return as wkt """
    ogrgeom = CreateGeometryFromWkt(shape)
    trans = CoordinateTransformation(SpatialReference(ssrs), SpatialReference(tsrs))
    ogrgeom.Transform(trans)
    wkt = ogrgeom.ExportToWkt()
    ogrgeom = None
    return wkt


def transform(filename, srs):
    """ Transform vector file to another SRS """
    # TODO - move functionality into GIPPY
    bname = os.path.splitext(os.path.basename(filename))[0]
    td = tempfile.mkdtemp()
    fout = os.path.join(td, bname + '_warped.shp')
    prjfile = os.path.join(td, bname + '.prj')
    f = open(prjfile, 'w')
    f.write(srs)
    f.close()
    cmd = 'ogr2ogr %s %s -t_srs %s' % (fout, filename, prjfile)
    result = commands.getstatusoutput(cmd)
    return fout


def crop2vector(img, vector):
    """ Crop a GeoImage down to a vector - only used by mosaic """
    # transform vector to srs of image
    vecname = transform(vector.Filename(), img.Projection())
    warped_vec = open_vector(vecname)
    # rasterize the vector
    td = tempfile.mkdtemp()
    mask = gippy.GeoImage(os.path.join(td, vector.LayerName()), img, gippy.GDT_Byte, 1)
    maskname = mask.Filename()
    mask = None
    cmd = 'gdal_rasterize -at -burn 1 -l %s %s %s' % (warped_vec.LayerName(), vecname, maskname)
    result = commands.getstatusoutput(cmd)
    VerboseOut('%s: %s' % (cmd, result), 4)
    mask = gippy.GeoImage(maskname)
    img.AddMask(mask[0]).Process().ClearMasks()
    mask = None
    shutil.rmtree(os.path.dirname(maskname))
    shutil.rmtree(os.path.dirname(vecname))
    # VerboseOut('Cropped to vector in %s' % (datetime.now() - start), 3)
    return img


def mosaic(images, outfile, vector):
    """ Mosaic multiple files together, but do not warp """
    nd = images[0][0].NoDataValue()
    srs = images[0].Projection()
    # check they all have same projection
    filenames = [images[0].Filename()]
    for f in range(1, images.NumImages()):
        if images[f].Projection() != srs:
            raise Exception("Input files have non-matching projections and must be warped")
        filenames.append(images[f].Filename())
    # transform vector to image projection
    geom = wktloads(transform_shape(vector.WKT(), vector.Projection(), srs))

    extent = geom.bounds
    ullr = "%f %f %f %f" % (extent[0], extent[3], extent[2], extent[1])

    # run merge command
    nodatastr = '-n %s -a_nodata %s -init %s' % (nd, nd, nd)
    cmd = 'gdal_merge.py -o %s -ul_lr %s %s %s' % (outfile, ullr, nodatastr, " ".join(filenames))
    result = commands.getstatusoutput(cmd)
    VerboseOut('%s: %s' % (cmd, result), 4)
    imgout = gippy.GeoImage(outfile, True)
    for b in range(0, images[0].NumBands()):
        imgout[b].CopyMeta(images[0][b])
    imgout.CopyColorTable(images[0])
    return crop2vector(imgout, vector)


def julian_date(date_and_time, variant=None):
    """Returns the julian date for the given datetime object.

    If no variant is chosen, the original julian date is given (days
    since noon, Jan 1, 4713 BC, fractions included).  If a variant is
    chosen, that julian date is returned instead.  Supported variants:
    'modified' (JD - 2400000.5) and 'cnes' (JD - 2433282.5).  See
    https://en.wikipedia.org/wiki/Julian_day for more details.
    """
    mjd_td = date_and_time - datetime.datetime(1858, 11, 17)
    # note day-length isn't constant under UTC due to leap seconds; hopefully this is close enough
    mjd = mjd_td.days + mjd_td.seconds / 86400.0

    offsets = {
        None:       2400000.5,
        'modified': 0.0,
        'cnes':     -33282.0,
    }

    return mjd + offsets[variant]


##############################################################################
# Error handling and script setup & teardown
##############################################################################

_traceback_verbosity = 4    # only print a traceback if the user selects this verbosity or higher
_accumulated_errors = []    # used for tracking success/failure & doing final error reporting when
                            # GIPS is running as a command-line application
_stop_on_error = False      # should GIPS try to recover from errors?  Set by gips_script_setup


def set_error_handler(handler):
    """Set the active error handler (generally for entire life of process)."""
    global error_handler
    error_handler = handler


def report_error(error, msg_prefix, show_tb=True):
    """Print an error report on stderr, possibly including a traceback.

    Caller can suppress the traceback with show_tb.  The user can suppress
    it via the GIPS global verbosity setting."""
    if show_tb and gippy.Options.Verbose() >= _traceback_verbosity:
        verbose_out(msg_prefix + ':', 1, stream=sys.stderr)
        traceback.print_exc()
    else:
        verbose_out(msg_prefix + ': ' + str(error), 1, stream=sys.stderr)


@contextmanager
def lib_error_handler(msg_prefix='Error', continuable=False):
    """Handle errors appropriately for GIPS running as a library."""
    try:
        yield
    except Exception as e:
        if continuable and not _stop_on_error:
            report_error(e, msg_prefix)
        else:
            raise


error_handler = lib_error_handler # set this so gips code can use the right error handler


def gips_exit():
    """Deliver an error report if needed, then exit."""
    if len(_accumulated_errors) == 0:
        sys.exit(0)
    verbose_out("Fatal: {} error(s) occurred:".format(len(_accumulated_errors)), 1, sys.stderr)
    [report_error(error, error.msg_prefix) for error in _accumulated_errors]
    sys.exit(1)


@contextmanager
def cli_error_handler(msg_prefix='Error', continuable=False):
    """Context manager for uniform error handling for command-line users.

    Exceptions are caught and reported to stderr; _gips_exit() is called
    if halt is indicated.
    """
    try:
        yield
    except Exception as e:
        e.msg_prefix = msg_prefix # for use by gips_exit
        _accumulated_errors.append(e)
        if continuable and not _stop_on_error:
            report_error(e, msg_prefix)
        else:
            gips_exit()


def gips_script_setup(driver_string=None, stop_on_error=False, setup_orm=True):
    """Run this at the beginning of a GIPS CLI program to do setup."""
    global _stop_on_error
    _stop_on_error = stop_on_error
    set_error_handler(cli_error_handler)
    from gips.inventory import orm # avoids a circular import
    with error_handler():
        # must run before orm.setup
        data_class = None if driver_string is None else import_data_class(driver_string)
        if setup_orm:
            orm.setup()
        return data_class
