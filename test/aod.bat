#!/bin/bash

ARGS="-s /etc/gips/test/NHseacoast.shp -d 2012-12-01,2012-12-2 -v 4"

gips_inventory AOD $ARGS --fetch
gips_process AOD $ARGS
gips_project AOD $ARGS --res 250 250 --outdir aod_project
gips_stats aod_project/0
gips_tiles AOD $ARGS --outdir aod_tiles
gips_stats aod_tiles

