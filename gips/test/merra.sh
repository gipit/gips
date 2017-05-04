#!/bin/bash

[ -z "$GIPSTESTPATH" ] && GIPSTESTPATH="."

ARGS="-s $GIPSTESTPATH/NHseacoast.shp -d 1980-01-01,1980-01-02 -v 4"

gips_info Merra 
gips_inventory merra $ARGS --fetch
gips_process merra $ARGS

# mosaic
gips_project merra $ARGS --res 100 100 --outdir merra_project --notld
gips_stats merra_project/*

# mosaic without warping
gips_project merra $ARGS --outdir merra_project_nowarp --notld
gips_stats merra_project_nowarp

# warp tiles
gips_tiles merra $ARGS --outdir merra_warped_tiles --notld
gips_stats merra_warped_tiles/*

# copy tiles
#gips_tiles merra -t h12v04 -d 2012-12-01,2012-12-10 -v 4 --outdir modis_tiles --notld
#gips_stats modis_tiles
