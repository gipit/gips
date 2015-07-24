#!/usr/bin/env python

try:
    from gips.data.cdl import cdlData as data
except:
    from gips.data.cdl import CDLData as data

from gips.core import SpatialExtent
from gips.utils import open_vector
from shapely.wkt import loads
from timeit import Timer

extent = SpatialExtent.factory(data, 'Belknap.shp')


assert len(extent[0].tiles) == 1, ('Belknap is contained wholly NH tile, '
                                   'but vector2tiles says it is in {} tiles'
                                   .format(len(extent[0].tiles)))

t = 'NH.shp'
f = 'NHseacoast.shp'


aWKT = open_vector(t)[0].WKT()
bWKT = open_vector(f)[0].WKT()
a = loads(aWKT)
b = loads(bWKT)
iters = 1000


def intersect():
    global a, b
    return a.intersects(b)


def overlap():
    global a, b
    return a.overlaps(b)


def relate():
    global a, b
    return a.relate(b).startswith('2')


def intersection():
    global a, b
    return a.intersection(b)


rel_total = Timer(relate).repeat(3, iters)
intersection_total = Timer(intersection).repeat(3, iters)
ovr_total = Timer(overlap).repeat(3, iters)
int_total = Timer(intersect).repeat(3, iters)


print('{} iterations\nintersects: {}\noverlaps: {}\nrelate.startwith("2"): {}\nintersection): {}'
      .format(iters, max(int_total), max(ovr_total), max(rel_total), max(intersection_total)))
