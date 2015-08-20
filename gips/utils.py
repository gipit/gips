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

import os
import errno
import gippy
from gippy import GeoVector
import tempfile
import commands
import shutil
import traceback


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


def VerboseOut(obj, level=1):
    if gippy.Options.Verbose() >= level:
        #pprint.PrettyPrinter().pprint(obj)
        if not isinstance(obj, (list, tuple)):
            obj = [obj]
        for o in obj:
            print o

##############################################################################
# Filesystem functions
##############################################################################

def File2List(filename):
    f = open(filename)
    txt = f.readlines()
    txt2 = []
    for t in txt:
        txt2.append(t.rstrip('\n'))
    return txt2


def List2File(lst, filename):
    f = open(filename, 'w')
    f.write('\n'.join(lst) + '\n')
    f.close()


def RemoveFiles(filenames, extensions=['']):
    for f in filenames:
        for ext in ([''] + extensions):
            try:
                os.remove(f + ext)
            except OSError as e:
                if e.errno != errno.ENOENT:
                    raise
                continue


def basename(str):
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

##############################################################################
# Settings functions
##############################################################################

def settings():
    """ Retrieve GIPS settings - first from user, then from system """
    import imp
    try:
        # import user settings first
        src = imp.load_source('settings', os.path.expanduser('~/.gips/settings.py'))
        return src
    except Exception, e:
        try:
            import gips.settings
            return gips.settings
        except:
            raise Exception('No settings found...did you run gips_config?')


def create_environment_settings(repos_path, email=''):
    """ Create settings file and data directory """
    from gips.settings_template import __file__ as src
    cfgpath = os.path.dirname(__file__)
    cfgfile = os.path.join(cfgpath, 'settings.py')
    try:
        if not os.path.exists(cfgfile):
            print 'Creating environment settings'
            with open(cfgfile, 'w') as fout:
                with open(src, 'r') as fin:
                    for line in fin:
                        fout.write(line.replace('$TLD', repos_path).replace('$EMAIL', email))
        return cfgfile
    except OSError:
        # no permissions, so no environment level config installed
        #print traceback.format_exc()
        return None


def create_user_settings(email=''):
    """ Create a settings file using the included template and the provided top level directory """
    from gips.user_settings_template import __file__ as src
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
    try:
        repos = settings().REPOS
    except:
        print traceback.format_exc()
        raise Exception('Problem reading repository...check settings files')
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
    found = False
    for key in sorted(repos.keys()):
        if os.path.isdir(repos[key]['repository']):
            try:
                repo = import_repository_class(key)
                sources[key] = repo.description
                found = True
            except:
                VerboseOut(traceback.format_exc(), 1)
        else:
            raise Exception('ERROR: archive %s is not a directory or is not available' % key)
    if not found:
        print "There are no available data sources!"
    return sources


def import_data_module(clsname):
    """ Import a data driver by name and return as module """
    import imp
    path = settings().REPOS[clsname].get('driver', '')
    if path == '':
        path = os.path.join( os.path.dirname(__file__), 'data', clsname, clsname + '.py') #__init__.py' )
    try:
        mod = imp.load_source(clsname, path)
        return mod
    except:
        print traceback.format_exc()


def import_repository_class(clsname):
    """ Get clsnameRepository class object """
    mod = import_data_module(clsname)
    exec('repo = mod.%sRepository' % clsname)
    return repo


def import_data_class(clsname):
    """ Get clsnameData class object """
    mod = import_data_module(clsname)
    exec('repo = mod.%sData' % clsname)
    return repo    


##############################################################################
# Geospatial functions
##############################################################################


from pdb import set_trace


def open_vector(fname, key="", where=''):
    """ Open vector or feature """
    parts = fname.split(':')
    if len(parts) == 1:
        vector = GeoVector(fname)
        vector.SetPrimaryKey(key)
    else:
        # or it is a database
        if parts[0] not in settings().DATABASES.keys():
            raise Exception("%s is not a valid database" % parts[0])
        try:
            db = settings().DATABASES[parts[0]]
            filename = ("PG:dbname=%s host=%s port=%s user=%s password=%s" %
                        (db['NAME'], db['HOST'], db['PORT'], db['USER'], db['PASSWORD']))
            vector = GeoVector(filename, parts[1])
            vector.SetPrimaryKey(key)

        except Exception, e:
            VerboseOut(traceback.format_exc(), 4)
    if where != '':
        # return array of features

        # return vector.where(where)
        # it looks like a GeoVector now only handles key=value and not
        # general SQL, which I think is less useful. So, temporary fix:
        key, value = where.split('=')
        return vector.where(key, value)

        features = []
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


# old code utilizing shared memory array
# Chunk it up
# chunksz = int(data.shape[0] / nproc)
# extra = data.shape[0] - chunksz * nproc
# chunks = [chunksz] * (nproc - extra) + [chunksz + 1] * extra

# queue = multiprocessing.Queue()
# from agspy.contrib import shmarray
# classmap = shmarray.create_copy(classmap)
# tmp = numpy.ctypeslib.as_ctypes(classmap)
# cmap = sharedctypes.Array(tmp._type_, tmp, lock=False)
