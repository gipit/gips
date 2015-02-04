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
from osgeo import ogr
from osgeo import osr
from shapely.wkb import loads
from shapely.geometry import Point
import datetime
import tempfile
import commands
#ogr.UseExceptions()

# OBSOLETE FILE - USE gippy.GeoVector 


# TODO - are both transform_point and transform_shape needed?
def transform_point(point, source, target):
    """ convert a single point (x, y) from one coordinate system to another """
    s_srs = osr.SpatialReference()
    s_srs.SetFromUserInput(source)
    t_srs = osr.SpatialReference()
    t_srs.SetFromUserInput(target)
    ct = osr.CoordinateTransformation(s_srs, t_srs)
    x, y = point
    newx, newy, z = ct.TransformPoint(x, y)
    return newx, newy


def transform_shape(shape, source, target):
    """ Convert a shapely shape to another SRS """
    srs_in = osr.SpatialReference()
    srs_out = osr.SpatialReference()
    srs_in.ImportFromWkt(source)
    srs_out.ImportFromWkt(target)
    trans = osr.CoordinateTransformation(srs_in, srs_out)
    ogrgeom = ogr.CreateGeometryFromWkb(shape.wkb)
    ogrgeom.Transform(trans)
    return loads(ogrgeom.ExportToWkb())


class GeoVector(object):
    """ A GeoVector object representing shapefile or PostGIS data layer """

    @property
    def name(self):
        return self.vector.GetName()

    @property
    def layer_name(self):
        return self.layer.GetName()

    @property
    def num_features(self):
        return self.layer.GetFeatureCount()

    @property
    def extent(self):
        ext = self.layer.GetExtent()
        return [ext[0], ext[2], ext[1], ext[3]]

    def proj(self):  # , format='Wkt'):
        format = 'Wkt'
        spatial_reference = self.layer.GetSpatialRef()
        return eval('spatial_reference.ExportTo%s()' % format)

    def identify(self, lon, lat):
        """ get information about a point """
        point = Point(lon, lat)
        for fid in self.get_fids():
            feature = self.get_geom(fid)
            if feature.contains(point):
                return self.get_feature(fid)
        return None

    def __init__(self, filename, layer=''):
        """ Open an existing vector file """
        self.filename = filename
        self.vector = ogr.Open(filename)
        if not self.vector:
            raise Exception("OGR can't open %s" % filename)

        # typical case is that shapefiles have one layer
        if layer == '':
            self.layer = self.vector.GetLayer(0)
        else:
            self.layer = self.vector.GetLayerByName(layer)

        # get list of attributes
        layer_defn = self.layer.GetLayerDefn()
        self.field_names = []
        for i in range(layer_defn.GetFieldCount()):
            self.field_names.append(layer_defn.GetFieldDefn(i).GetName())

    def union(self):
        """ Compute the union of all geometries in layer and return Shapely object """
        from shapely.ops import unary_union
        shapes = []
        #for i in range(self.num_features):
        #    feat = self.layer.GetFeature(i)
        for feature in self.layer:
            shapes.append(loads(feature.GetGeometryRef().ExportToWkb()))
        return unary_union(shapes)

    def transform(self, srs, filename=''):
        """ Transform to another SRS and return """
        bname = os.path.splitext(os.path.basename(self.filename))[0]
        td = tempfile.mkdtemp()
        if filename == '':
            filename = os.path.join(td, bname+'_warped.shp')
        prjfile = os.path.join(td, bname+'.prj')
        f = open(prjfile, 'w')
        f.write(srs)
        f.close()
        cmd = 'ogr2ogr %s %s -t_srs %s' % (filename, self.filename, prjfile)
        result = commands.getstatusoutput(cmd)
        return GeoVector(filename)

    def get_fids(self):
        """ Return list of feature IDs """
        feature_list = []
        for ifeature in range(self.num_features):
            feature = self.layer.GetFeature(ifeature)
            fid = feature.GetFID()
            assert fid == ifeature
            feature_list.append(fid)
        return feature_list

    def get_geom(self, fid):
        """ Get shapely geometry for specified fid """
        feature = self.layer.GetFeature(fid)
        return loads(feature.GetGeometryRef().ExportToWkb())

    def get_feature(self, fid):
        """ get all the attributes for a specified feature (get a row) """
        feature = self.layer.GetFeature(fid)
        result = {}
        for fieldname in self.field_names:
            field_index = feature.GetFieldIndex(fieldname)
            result[fieldname] = feature.GetField(field_index)
        return result

    def get_attributes(self, fieldname):
        """ get an attribute for all features (get a column) """
        column = []
        for ifeature in range(self.num_features):
            feature = self.layer.GetFeature(ifeature)
            fid = feature.GetFID()
            assert fid == ifeature
            field_index = feature.GetFieldIndex(fieldname)
            column.append(feature.GetField(field_index))
        return column

    def get_attribute(self, fid, fieldname):
        """ get an attribute for a specified feature and column """
        feature = self.layer.GetFeature(fid)
        field_index = feature.GetFieldIndex(fieldname)
        result = feature.GetField(field_index)
        return result

    def get_data(self, as_numpy=False):
        """ get all data as a table """
        assert '_x' not in self.field_names
        assert '_y' not in self.field_names
        result = {'_x': [], '_y': []}
        for field_name in self.field_names:
            result[field_name] = []
        for ifeature in range(self.num_features):
            feature = self.layer.GetFeature(ifeature)
            fid = feature.GetFID()
            assert fid == ifeature
            try:
                # point feature
                geom = feature.GetGeometryRef()
                point = geom.GetPoint()[:2]
                result['_x'].append(point[0])
                result['_y'].append(point[1])
            except:
                # polygon feature
                #point = self.feature_centroid(fid)
                #print 'polygon'
                points = self.feature_vertices(fid)
                xpts = [point[0] for point in points]
                ypts = [point[1] for point in points]
                result['_x'].append(xpts)
                result['_y'].append(ypts)
            for field_name in self.field_names:
                field_index = feature.GetFieldIndex(field_name)
                result[field_name].append(feature.GetField(field_index))
        if as_numpy:
            import numpy as np
            for k, v in result.items():
                result[k] = np.array(v)
        return result

    # TODO - validate these next two functions
    def feature_extent(self, fid, srs=None, pad=0):
        """ calculate geographic extent of a feature """
        feature = self.layer.GetFeature(fid)
        geom = feature.GetGeometryRef()
        ring = geom.GetGeometryRef(0)
        npoints = ring.GetPointCount()
        xmin, xmax, ymin, ymax = [9.e99, -9.e99, 9.e99, -9.e99]
        for ipoint in xrange(npoints):
            xpt, ypt, zpt = ring.GetPoint(ipoint)
            if srs:
                xpt, ypt = transform_point((xpt, ypt), self.proj_wkt, srs)
            if xpt < xmin:
                xmin = xpt
            if ypt < ymin:
                ymin = ypt
            if xpt > xmax:
                xmax = xpt
            if ypt > ymax:
                ymax = ypt
        return xmin, ymin, xmax, ymax

    def feature_centroid(self, fid, srs=None):
        """ calculate centroid of a feature """
        feature = self.layer.GetFeature(fid)
        geom = feature.GetGeometryRef()
        ring = geom.GetGeometryRef(0)
        npoints = ring.GetPointCount()
        xmin, xmax, ymin, ymax = [9.e99, -9.e99, 9.e99, -9.e99]
        for ipoint in xrange(npoints):
            xpt, ypt, zpt = ring.GetPoint(ipoint)
            if xpt < xmin:
                xmin = xpt
            if ypt < ymin:
                ymin = ypt
            if xpt > xmax:
                xmax = xpt
            if ypt > ymax:
                ymax = ypt
        xpt = (xmax + xmin)/2.
        ypt = (ymax + ymin)/2.
        if srs:
            xpt, ypt = transform_point((xpt, ypt), self.proj_wkt, srs)
        return xpt, ypt

    def polycontainspoints(self, shapefile, pointfile):
        """
        Take comma-delimited lon,lat points from a pointfile and test each one for membership
        within the polygon specified by the shapefile.
        """
        from shapely.wkb import loads
        from shapely.geometry import Point
        # Open the shapefile
        source = ogr.Open(shapefile)
        # Extract the first layer, assume it is the only one
        layer = source.GetLayer(0)
        # Get the first feature, assume it is the only one
        feature = layer.GetNextFeature()
        # Convert the OGR polygon into a Shapely polygon using WKB (Well-Known Binary) format
        polygon = loads(feature.GetGeometryRef().ExportToWkb())
        # Read the lon,lat points from the file
        lonlats = open(pointfile, "r").readlines()
        # Initialize the result array
        result = []
        # Loop over the points, there's a faster way to do this, see Shapely manual section 5.1.1
        for lonlat in lonlats:
            lonlat = lonlat.split(",")
            lon, lat = [float(ll) for ll in lonlat]
            point = Point(lon, lat)
            within = polygon.contains(point)
            result.append((lon, lat, within))
        # Give back the result
        return result

    # TODO - refactor below methods. clean up data type functions
    def add_attribute(self, val_column, col_name, datatype=str):
        """ Add colum data to attribute table """
        # TODO: FIX
        dtype = get_column_type(val_column[col_name])
        ogr_dtype = get_ogr_type(dtype)
        self.layer.CreateField(ogr.FieldDefn(col_name, ogr_dtype))
        for ifeature in range(self.num_features):
            feature = self.layer.GetFeature(ifeature)
            self.layer.DeleteFeature(ifeature)
            feature.SetField(col_name, val_column[ifeature])
            self.layer.CreateFeature(feature)
        return None

    def join(self, val_columns, joinkey, outpath):
        """ create additional attributes based on external data """
        joincol = val_columns.pop(joinkey)
        colnames = val_columns.keys()
        newshp = NewShapefile(outpath, layer=self.layer)
        for colname in colnames:
            dtype = get_column_type(val_columns[colname])
            ogr_dtype = get_ogr_type(dtype)
            newshp.layer.CreateField(ogr.FieldDefn(colname, ogr_dtype))
        nfeatures = newshp.layer.GetFeatureCount()
        for ifeature in range(nfeatures):
            feature = newshp.layer.GetFeature(ifeature)
            joinloc = feature.GetFieldIndex(joinkey)
            joinval = feature.GetField(joinloc)
            if joinval is None:
                continue
            joinidx = joincol.index(joinval)
            # note val_columns is a dictionary
            for colname, val_column in val_columns.items():
                feature.SetField(colname, val_column[joinidx])
            newshp.layer.DeleteFeature(ifeature)
            newshp.layer.CreateFeature(feature)
            feature = None
        return None

    # TODO - refactor below classmethods, provide factory functions for points, lines, polys
    @classmethod
    def NewShapeFile(cls, filename, srs=None):
        """ Create new polygon shapefile """
        drv = ogr.GetDriverByName('ESRI Shapefile')
        if os.path.exists(filename):
            os.remove(filename)
        vector = drv.CreateDataSource(filename)
        proj = osr.SpatialReference()
        proj.SetFromUserInput(srs)
        layer = vector.CreateLayer(vector.GetName(), geom_type=ogr.wkbPolygon, srs=self.proj)
        vector = None
        return GeoVector(filename)

    @classmethod
    def NewShapeFile_point(cls, outpath, data, proj='EPSG:4326', lonkey='x', latkey='y', maxchar=80):
        """ Create new point shapefile """
        if os.path.exists(outpath):
            os.remove(outpath)
        drv = ogr.GetDriverByName('ESRI Shapefile')
        shape = drv.CreateDataSource(outpath)
        srs = osr.SpatialReference()
        srs.SetFromUserInput(proj)
        name = os.path.splitext(os.path.split(outpath)[1])[0]
        layer = shape.CreateLayer(name, geom_type=ogr.wkbPoint, srs=srs)
        colnames = data.keys()
        keys = data.keys()
        for i, colname in enumerate(colnames):
            colnames[i] = colname[:10]
        nfeatures = len(data[keys[0]])
        for i, key in enumerate(keys):
            if key in (lonkey, latkey):
                continue
            dtype = get_column_type(data[key])
            ogr_dtype = get_ogr_type(dtype)
            field_def = ogr.FieldDefn(colnames[i], ogr_dtype)
            if ogr_dtype == ogr.OFTString:
                field_def.SetWidth(maxchar)
            layer.CreateField(field_def)
        feature = ogr.Feature(feature_def=layer.GetLayerDefn())
        geom = ogr.Geometry(type=ogr.wkbPoint)
        for ifeature in range(nfeatures):
            try:
                x = data[lonkey][ifeature]
                y = data[latkey][ifeature]
                geom.SetPoint(0, float(x), float(y))
            except:
                print "skipping feature", ifeature, x, y
                continue
            feature.SetGeometry(geom)
            feature.SetFID(ifeature)
            for i, key in enumerate(keys):
                if key in (lonkey, latkey):
                    continue
                feature.SetField(colnames[i], data[key][ifeature])
            layer.CreateFeature(feature)


def get_column_type(col):
    """ detect best python type to use for column """
    for i, item in enumerate(col):
        item = retype(item)
        if i == 0:
            coltype = type(item)
        newtype = type(item)
        if coltype is not newtype:
            coltype = getrighttype(newtype, coltype)
    return coltype


def get_ogr_type(value):
    """ map OGR types to python types """
    lookup = {str: ogr.OFTString, float: ogr.OFTReal, int: ogr.OFTInteger,
              datetime: ogr.OFTDateTime, date: ogr.OFTDate, time: ogr.OFTTime}
    # 4, 2, 0, 11, 9, 10
    if type(value) is type:
        return lookup[value]
    else:
        return lookup[type(value)]


def getrighttype(newtype, oldtype):
    """ define compatible types """
    if newtype is str or oldtype is str:
        return str
    if newtype is int and oldtype is float:
        return float
    if newtype is float and oldtype is int:
        return float
    if newtype is datetime.datetime and oldtype is not datetime.datetime:
        return str
    if newtype is not datetime.datetime and oldtype is datetime.datetime:
        return str
    if newtype is datetime.date and oldtype is not datetime.date:
        return str
    if newtype is not datetime.date and oldtype is datetime.date:
        return str
    if newtype is datetime.time and oldtype is not datetime.time:
        return str
    if newtype is not datetime.time and oldtype is datetime.time:
        return str


def retype(arg, datestr='%Y-%m-%d'):
    """ attempt to recast a string as int, float, or datetime """
    charstr = str(arg).strip()
    try:
        value = int(arg)
    except:
        try:
            value = float(arg)
        except:
            try:
                value = datetime.datetime.strptime(arg, datestr)
            except:
                value = arg
    return value

"""
def intersection(GeoVector1, GeoVector2):
    geom1 = GeoVector1.union()
    ogrgeom = ogr.CreateGeometryFromWkb(geom1.wkb)
    GeoVector2.layer.SetSpatialFilter(ogrgeom)
    GeoVector2.layer.ResetReading()
    feat = GeoVector2.layer.GetNextFeature()
    fldindex = feat.GetFieldIndex('pr')
    dict = {}  # 'total':geom1.area}
    while feat is not None:
        geom2 = loads(feat.GetGeometryRef().ExportToWkb())
        area = geom1.intersection(geom2).area
        if area != 0:
            dict[feat.GetField(fldindex)] = area/geom1.area
        feat = GeoVector2.layer.GetNextFeature()
    return dict
    #tmpdir = tempfile.mkdtemp()
    #vout = NewShapefile("%s/intersection.shp" % tmpdir, srs=GeoVector1.get_proj())
    #GeoVector1.layer.Intersection(GeoVector2.layer, vout.layer)
    #return vout
"""