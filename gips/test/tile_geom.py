#!/usr/bin/env python

from gips.data.cdl import cdlRepository as repo, cdlData as data
from gips.core import SpatialExtent

extent = SpatialExtent.factory(data, 'Belknap.shp')


assert len(extent[0].tiles) == 1, ('Belknap is contained wholly NH tile, '
                                   'but vector2tiles says it is in {} tiles'
                                   .format(len(extent[0].tiles)))
