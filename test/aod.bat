#!/bin/bash

ARGS="-s /etc/gips/test/NHseacoast.shp -d 2012-12-01,2012-12-2 -v 4"

gips_info AOD
gips_inventory AOD $ARGS --fetch

gips_process AOD $ARGS

# mosaic
gips_project AOD $ARGS --res 250 250 --outdir aod_project --notld
gips_stats aod_project/*

# warp tiles
gips_tiles AOD $ARGS --outdir aod_warped_tiles --notld
gips_stats aod_warped_tiles/*

# copy tiles
gips_tiles AOD -d 2012-12-01,2012-12-10 --outdir aod_tiles --notld
gips_stats aod_tiles
