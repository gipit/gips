#!/bin/bash

# Modis tests
ARGS="-s /etc/gips/test/NH.shp -d 2012-12 -v 3"
gips_inventory Modis $ARGS --fetch
gips_process Modis $ARGS -p ndvi8 indices
