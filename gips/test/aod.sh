#!/bin/bash

set -e

[ -z "$GIPSTESTPATH" ] && GIPSTESTPATH="."

ARGS="-s $GIPSTESTPATH/NHseacoast.shp -d 2012-12-01,2012-12-2 -v 4"

gips_info aod
gips_inventory aod $ARGS --fetch

gips_process aod $ARGS

# mosaic
gips_project aod $ARGS --res 250 250 --outdir aod_project --notld
gips_stats aod_project/*

# mosaic without warping
gips_project aod $ARGS --outdir aod_project_nowarp --notld
gips_stats aod_project_nowarp

# warp tiles
gips_tiles aod $ARGS --outdir aod_warped_tiles --notld
gips_stats aod_warped_tiles/*

# copy tiles
gips_tiles aod -d 2012-12-01,2012-12-10 --outdir aod_tiles --notld
gips_stats aod_tiles
