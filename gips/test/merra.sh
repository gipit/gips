#!/bin/bash

ARGS="-s /etc/gips/test/NHseacoast.shp -d 2012-12-01,2012-12-10 -v 4"

gips_info Merra 
gips_inventory Merra $ARGS --fetch
gips_process Merra $ARGS

# mosaic
gips_project Merra $ARGS --res 100 100 --outdir merra_project --notld
gips_stats merra_project/*

# mosaic without warping
gips_project Merra $ARGS --outdir merra_project_nowarp --notld
gips_stats merra_project_nowarp

# warp tiles
gips_tiles Merra $ARGS --outdir merra_warped_tiles --notld
gips_stats merra_warped_tiles/*

# copy tiles
#gips_tiles Merra -t h12v04 -d 2012-12-01,2012-12-10 -v 4 --outdir modis_tiles --notld
#gips_stats modis_tiles
