#!/bin/bash

ARGS="-s /etc/gips/test/NHseacoast.shp -d 2012-12 -v 3 -p ndvi8"

# Modis tests
PRODUCTS="-p ndvi8"
gips_info Modis
gips_inventory Modis $ARGS --fetch
gips_process Modis $ARGS
gips_project Modis $ARGS --res 100 100 --outdir modis_project
gips_stats modis_project/0
gips_tiles Modis $ARGS --outdir modis_tiles
gips_stats modis_tiles

# AOD tests
#gips_inventory AOD $ARGS --fetch
#gips_process AOD $ARGS
#gips_project AOD $ARGS --res 250 250 --outdir aod_project
#gips_stats aod_project/0
#gips_tiles Modis $ARGS --outdir aod_tiles
#gips_stats aod_tiles

