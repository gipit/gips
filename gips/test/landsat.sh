#!/bin/bash

set -e

DATES="-d 2012-256"

[ -z "$GIPSTESTPATH" ] && GIPSTESTPATH="."

ARGS="-s $GIPSTESTPATH/NHseacoast.shp $DATES -v 4 -p acca ref-toa ndvi-toa rad-toa"

gips_inventory landsat $ARGS
gips_process landsat $ARGS

# mosaic
gips_project landsat $ARGS --res 30 30 --outdir landsat_project --notld
gips_stats landsat_project/*

# mosaic without warping
gips_project landsat $ARGS --outdir landsat_project_nowarp --notld
gips_stats landsat_project_nowarp

# warp tiles
gips_tiles landsat $ARGS --outdir landsat_warped_tiles --notld
gips_stats landsat_warped_tiles/*

# copy tiles
gips_tiles landsat -t 012030 $DATES -v 4 --outdir landsat_tiles --notld -p ref-toa ndvi-toa rad-toa
gips_stats landsat_tiles
