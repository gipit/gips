#!/usr/bin/env python

import os
import glob
from osgeo import gdal

import numpy as np
import fiona
from fiona.crs import from_string
from shapely.geometry import mapping, Polygon


from pdb import set_trace


STARTDIR = "/titan/data/modis/tiles"
DATE = '2004001'
WIDTH = 2400
OUTPATH = "new/tiles.shp"
PROJ = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"


def write_shapefile(data, proj):
    schema = {
        'geometry': 'Polygon',
        'properties': {'tileid': 'str', 'bounds': 'str'},
    }
    crs = from_string(proj)
    with fiona.open(OUTPATH, 'w', 'ESRI Shapefile', schema, crs=crs) as shp:
        for tileid, coords in data:

            arr = np.array(coords)
            xmin, ymin = tuple(arr.min(axis=0))
            xmax, ymax = tuple(arr.max(axis=0))
            bounds = str((xmin, ymin, xmax, ymax))

            poly = Polygon(coords)
            shp.write({
                'geometry': mapping(poly),
                'properties': {'tileid': tileid, 'bounds': bounds},
            })


def main():
    tiles = glob.glob(os.path.join(STARTDIR, '*'))
    itile = 0
    boxes = []
    for tile in tiles:
        lcglob = os.path.join(tile, DATE, '*landcover*')
        try:
            lcpath = glob.glob(lcglob)[0]
            print lcpath
        except:
            continue
        gdalfile = gdal.Open(lcpath)
        name = os.path.splitext(os.path.split(lcpath)[1])[0]
        tile = name.split('_')[0]
        x0, dx, dxy, y0, dyx, dy = gdalfile.GetGeoTransform()
        xs, ys = (x0, y0)
        coords = [(x0, y0)]
        for direction in ((1, 0), (0, 1), (-1, 0), (0, -1)):
            for step in range(10):
                x1 = x0 + 240*dx*direction[0]
                y1 = y0 + 240*dy*direction[1]
                coords.append((x1, y1))
                x0, y0 = (x1, y1)
        coords[-1] = (xs, ys)
        print tile, coords
        boxes.append((tile, coords))
        itile += 1

    write_shapefile(boxes, PROJ)


if __name__ == "__main__":
    main()
