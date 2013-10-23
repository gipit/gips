#!/usr/bin/env python

import sys
import ogr, osr
from shapely.wkb import loads
from shapely.geometry import Point
from datetime import datetime, date, time
import tempfile

ogr.UseExceptions()

def intersection(GeoVector1, GeoVector2):
    geom1 = GeoVector1.union()
    ogrgeom = ogr.CreateGeometryFromWkb(geom1.wkb)
    print "Site Area = ",geom1.area
    GeoVector2.layer.SetSpatialFilter(ogrgeom)
    GeoVector2.layer.ResetReading()
    feat = GeoVector2.layer.GetNextFeature()
    fldindex = feat.GetFieldIndex('pr')
    dict = {} #'total':geom1.area}
    while feat is not None:
        geom2 = loads(feat.GetGeometryRef().ExportToWkb())
        area = geom1.intersection(geom2).area
        if area != 0: dict[feat.GetField(fldindex)] = area/geom1.area
        feat = GeoVector2.layer.GetNextFeature()
    return dict
    #tmpdir = tempfile.mkdtemp()
    #vout = NewShapefile("%s/intersection.shp" % tmpdir, srs=GeoVector1.get_proj())
    #GeoVector1.layer.Intersection(GeoVector2.layer, vout.layer)
    #return vout

def transform_point(point, source="WGS84", target="WGS84"):
    """ convert a single point (x, y) from one coordinate system to another """
    s_srs = osr.SpatialReference()
    s_srs.SetFromUserInput(source)
    t_srs = osr.SpatialReference()
    t_srs.SetFromUserInput(target)    
    ct = osr.CoordinateTransformation(s_srs, t_srs)
    x, y = point
    newx, newy, z = ct.TransformPoint(x, y)    
    return newx, newy

def get_ogr_type(value):
    """ map OGR types to python types """
    lookup = {str: ogr.OFTString, float: ogr.OFTReal, int: ogr.OFTInteger, 
                datetime: ogr.OFTDateTime, date: ogr.OFTDate, time: ogr.OFTTime}
    # 4, 2, 0, 11, 9, 10
    if type(value) is type:
        return lookup[value]
    else:
        return lookup[type(value)]
  
def get_column_type(col):
    """ detect best python type to use for column """
    for i,item in enumerate(col):
        item = retype(item)
        if i==0:
            coltype = type(item)
        newtype = type(item)
        if coltype is not newtype:
            coltype = getrighttype(newtype, coltype)        
    return coltype

def getrighttype(newtype, oldtype):
    """ define compatible types """
    if newtype is str or oldtype is str:
        return str
    if newtype is int and oldtype is float:
        return float
    if newtype is float and oldtype is int:
        return float
    if newtype is datetime and oldtype is not datetime:
        return str
    if newtype is not datetime and oldtype is datetime:
        return str
    if newtype is date and oldtype is not date:
        return str
    if newtype is not date and oldtype is date:
        return str
    if newtype is time and oldtype is not time:
        return str
    if newtype is not time and oldtype is time:
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
                value = datetime.strptime(arg, datestr)
            except:
                value = arg
    return value


def polycontainspoints(shapefile, pointfile):
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
    lonlats = open(pointfile,"r").readlines()
    # Initialize the result array
    result = []
    # Loop over the points, there's a faster way to do this, see Shapely manual section 5.1.1
    for lonlat in lonlats:
        lonlat = lonlat.split(",")
        lon,lat = [float(ll) for ll in lonlat]
        point = Point(lon, lat)
        within = polygon.contains(point)
        result.append((lon, lat, within))
    # Give back the result
    return result


class GeoVector(object):
    def __init__(self, filename, layer=''):
        self.shape = ogr.Open(filename)
        if not self.shape:
            raise Exception, "OGR can't open %s" % filename
        self.nlayers = self.shape.GetLayerCount()
        self.name = self.shape.GetName()
        # typical case is that shapefiles have one layer
        if layer == '':
            self.layer = self.shape.GetLayer(0)
        else:
            self.layer = self.shape.GetLayerByName(layer)
        self.layername = self.layer.GetName()
        self.nfeatures = self.layer.GetFeatureCount()
        self.proj_wkt = self.get_proj('Wkt')
        layer_defn = self.layer.GetLayerDefn()
        self.field_names = []
        for i in range(layer_defn.GetFieldCount()):
            self.field_names.append(layer_defn.GetFieldDefn(i).GetName())
        self.extent = self.layer.GetExtent()

    def get_proj(self, format='Wkt'):
        spatial_reference = self.layer.GetSpatialRef()
        export_function = eval('spatial_reference.ExportTo'+format)
        return export_function()

    def copy(self, new_name='new_name'):
        """ make a copy of the shapefile """
        return NewShapefile(new_name, layer=self.layer)

    def add_column(self, val_column, col_name, datatype=str):
        """ append data to attribute table """
        # TODO: FIX
        dtype = get_column_type(val_column[col_name])
        ogr_dtype = get_ogr_type(dtype)
        self.layer.CreateField(ogr.FieldDefn(col_name, ogr_dtype))
        for ifeature in range(self.nfeatures):
            feature = self.layer.GetFeature(ifeature)
            self.layer.DeleteFeature(ifeature)
            feature.SetField(col_name, val_column[ifeature])
            self.layer.CreateFeature(feature)
        return None

    def union(self):
        """ Compute the union of all geometries in layer """
        from shapely.ops import unary_union
        shapes = []
        for i in range(self.nfeatures):
            feat = self.layer.GetFeature(i)
            shapes.append(loads(feat.GetGeometryRef().ExportToWkb()))
        return unary_union(shapes)

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
            if joinval is None: continue
            joinidx = joincol.index(joinval)
            # note val_columns is a dictionary
            for colname, val_column in val_columns.items():
                feature.SetField(colname, val_column[joinidx])
            newshp.layer.DeleteFeature(ifeature)
            newshp.layer.CreateFeature(feature)
            feature = None
        return None
        
    def feature_attributes(self, fid):
        """ get all the attributes for a specified feature (get a row) """
        feature = self.layer.GetFeature(fid)
        result = {}
        for fieldname in self.field_names:
            field_index = feature.GetFieldIndex(fieldname)
            result[fieldname] = feature.GetField(field_index)
        return result

    def get_geom(self, fid):
        """ Get shapely geometry for specified fid """
        feature = self.layer.GetFeature(fid)
        return loads(feature.GetGeometryRef().ExportToWkb())

    def get_attribute(self, fid, fieldname):
        """ get an attribute for a specified feature and column """
        feature = self.layer.GetFeature(fid)
        field_index = feature.GetFieldIndex(fieldname)
        result = feature.GetField(field_index)
        return result

    def attribute_column(self, field_name):
        """ get an attribute for all features (get a column) """
        column = []
        for ifeature in range(self.nfeatures):
            feature = self.layer.GetFeature(ifeature)
            fid = feature.GetFID()
            assert fid == ifeature
            field_index = feature.GetFieldIndex(field_name)
            column.append(feature.GetField(field_index))
        return column

    def extract_data(self, as_numpy=False):
        """ get all data as a table """
        assert '_x' not in self.field_names
        assert '_y' not in self.field_names
        result = {'_x':[], '_y':[]}
        for field_name in self.field_names:
            result[field_name] = []
        for ifeature in range(self.nfeatures):
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
            for k,v in result.items():
                result[k] = np.array(v)
        return result

    def feature_extent(self, fid, t_srs=None, pad=0, order='xyxy'):
        """ calculate geographic extent of a feature """
        assert order in ('xxyy', 'xyxy')
        feature = self.layer.GetFeature(fid)
        geom = feature.GetGeometryRef()
        ring = geom.GetGeometryRef(0)
        npoints = ring.GetPointCount()
        xmin, xmax, ymin, ymax = [9.e99, -9.e99, 9.e99, -9.e99]
        for ipoint in xrange(npoints):
            xpt, ypt, zpt = ring.GetPoint(ipoint)
            if t_srs:
                xpt, ypt = transform_point((xpt, ypt), self.proj_wkt, t_srs)
            if xpt < xmin: xmin = xpt
            if ypt < ymin: ymin = ypt
            if xpt > xmax: xmax = xpt
            if ypt > ymax: ymax = ypt
        if pad:
            xmin -= pad*(xmax - xmin)
            xmax += pad*(xmax - xmin)
            ymin -= pad*(ymax - ymin)
            ymax += pad*(ymax - ymin)
        if order == 'xyxy':
            return xmin, ymin, xmax, ymax
        else:
            return xmin, xmax, ymin, ymax

    def feature_centroid(self, fid, t_srs=None):
        """ calculate centroid of a feature """
        feature = self.layer.GetFeature(fid)
        geom = feature.GetGeometryRef()
        ring = geom.GetGeometryRef(0)
        npoints = ring.GetPointCount()
        xmin, xmax, ymin, ymax = [9.e99, -9.e99, 9.e99, -9.e99]
        for ipoint in xrange(npoints):
            xpt, ypt, zpt = ring.GetPoint(ipoint)
            if xpt < xmin: xmin = xpt
            if ypt < ymin: ymin = ypt
            if xpt > xmax: xmax = xpt
            if ypt > ymax: ymax = ypt
        xpt = (xmax + xmin)/2.
        ypt = (ymax + ymin)/2.        
        if t_srs:
            xpt, ypt = transform_point((xpt, ypt), self.proj_wkt, t_srs)
        return xpt, ypt

    def feature_vertices(self, fid, t_srs=None):
        feature = self.layer.GetFeature(fid)
    
        geom = feature.GetGeometryRef()        
        npoints = geom.GetGeometryCount()

        for ipoint in xrange(npoints):
            ring = geom.GetGeometryRef(ipoint)
            #print ring.ExportToJson()
            coords = eval(ring.ExportToJson())['coordinates'][0]
            print len(coords)
                        
        sys.exit()

        """        
        pts = []
        for ipoint in xrange(npoints):
            xpt, ypt, zpt = ring.GetPoint(ipoint)
            if t_srs:
                xpt, ypt = transform_point((xpt, ypt), self.proj_wkt, t_srs)
            pts.append((xpt, ypt))
        return pts
        """

    def identify(self, lon, lat):
        """ get information about a point """
        point = Point(lon, lat)
        feat_info = None
        for ifeature in range(self.nfeatures):
            feature = self.layer.GetFeature(ifeature)
            fid = feature.GetFID()
            assert fid == ifeature
            # convert the OGR polygon into a Shapely polygon using WKB
            polygon = loads(feature.GetGeometryRef().ExportToWkb())
            if polygon.contains(point):
                feat_info = self.feature_attributes(fid)
                break
        return feat_info

    def get_fids(self):
        """ extract IDs from shape file """
        feature_list = []
        for ifeature in range(self.nfeatures):
            feature = self.layer.GetFeature(ifeature)
            fid = feature.GetFID()
            assert fid == ifeature
            feature_list.append(fid)
        self.layer = self.shape.GetLayer(0)
        return feature_list

    def features(self):
        for ifeature in range(self.nfeatures):
            feature = self.layer.GetFeature(ifeature)
            fid = feature.GetFID()
            assert fid == ifeature
            yield fid, feature

class NewShapefile(object):
    def __init__(self, outpath, srs=None, layer=None, data=None):
        import os
        # only specify one of layer or srs
        assert bool(layer) ^ bool(srs)
        drv = ogr.GetDriverByName('ESRI Shapefile')
        if os.path.exists(outpath):
            os.remove(outpath)
        self.shape = drv.CreateDataSource(outpath)
        self.name = self.shape.GetName()
        if not layer:
            self.proj = osr.SpatialReference()
            self.proj.SetFromUserInput(srs)
            self.layer = self.shape.CreateLayer(self.name, geom_type=ogr.wkbPolygon, srs=self.proj)
        else:
            self.shape.CopyLayer(layer, self.name)
            self.layer = self.shape.GetLayer(0)

def write_pt_shpfile(outpath, data, proj='EPSG:4326', lonkey='x', latkey='y', maxchar=80):
    """ create and output a new point shapefile """
    import os
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
    for i,colname in enumerate(colnames):
        colnames[i] = colname[:10]
    nfeatures = len(data[keys[0]])
    for i,key in enumerate(keys):
        if key in (lonkey, latkey): continue
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
        for i,key in enumerate(keys):
            if key in (lonkey, latkey): continue
            feature.SetField(colnames[i], data[key][ifeature])
        layer.CreateFeature(feature)
   

        
        
 
