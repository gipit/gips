#!/bin/bash

ARGS="-s /etc/gips/test/NHseacoast.shp -d 2012-12-01,2012-12-10 -v 4"

gips_info Modis
gips_inventory Modis $ARGS --fetch
gips_process Modis $ARGS
gips_project Modis $ARGS --res 100 100 --outdir modis_project --notld
gips_stats modis_project/*

# warp tiles
gips_tiles Modis $ARGS --outdir modis_warped_tiles --notld
gips_stats modis_warped_tiles/*

# copy tiles
gips_tiles Modis -t h12v04 -d 2012-12-01,2012-12-10 -v 4 --outdir modis_tiles --notld
gips_stats modis_tiles
