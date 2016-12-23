
from osgeo import gdal, osr, ogr


def write_raster(fname, data, proj, geo, meta, bandnames=[], gcps=None, nodata=None):
    driver = gdal.GetDriverByName('GTiff')
    try:
        (nband, ny, nx) = data.shape
    except:
        # TODO error-handling-fix: leave as-is but report the error
        nband = 1
        (ny, nx) = data.shape
        data = data.reshape(1, ny, nx)
    dtype = create_datatype(data.dtype)
    tfh = driver.Create(fname, nx, ny, nband, dtype, [])
    tfh.SetGeoTransform(geo)
    tfh.SetMetadata(meta)
    if gcps is None:
        tfh.SetProjection(proj)
    else:
        tfh.SetGCPs(gcps, proj)
    if bandnames:
        assert len(bandnames) == nband
    for i in range(nband):
        band = tfh.GetRasterBand(i+1)
        if bandnames:
            assert len(bandnames) == nband, "wrong number of band names"
            print i, bandnames[i]
            band.SetDescription(bandnames[i])
        if nodata is not None:
            band.SetNoDataValue(nodata)
        band.WriteArray(data[i])
    tfh = None

gdal.GDT_UInt8 = gdal.GDT_Byte

def create_datatype(np_dtype):
    """ provide translation between data type codes """
    np_dtype = str(np_dtype)
    typestr = 'gdal.GDT_' + np_dtype.title().replace('Ui', 'UI')
    g_dtype = eval(typestr)
    return g_dtype

def create_meta(descr):
    """ generate metadata based on a simple template """
    timestamp = datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
    meta = {'AGS_Description':descr, 'AGS_Created':timestamp}
    return meta

def create_proj(crs, source='EPSG'):
    """ provide WKT formatted projection information based on an EPSG code """
    ref = osr.SpatialReference()
    eval('ref.ImportFrom%s(crs)' % source)
    return ref.ExportToWkt()

def create_geo(xycorner, resx, resy=None, northup=True):
    """ generate spatial reference information """
    x0, y0 = xycorner
    if resy is None:
        dx, dy = (resx, resx)
    else:
        dx, dy = (resx, resy)
    dxdy, dydx = (0, 0)
    if northup and resy is None: dy = -dy
    geo = (x0, dx, dxdy, y0, dydx, dy)
    return geo

