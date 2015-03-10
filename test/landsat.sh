#!/bin/bash

SHAPE="-s /etc/gips/test/NHseacoast.shp"
DATES="-d 2012-256"

ARGS="-s /etc/gips/test/NHseacoast.shp -d 2012-256 -v 4 -p ref-toa ndvi-toa rad-toa"

gips_info Landsat
gips_inventory Landsat $ARGS
gips_process Landsat $ARGS

# mosaic
gips_project Landsat $ARGS --res 30 30 --outdir landsat_project --notld
gips_stats landsat_project/*

# mosaic without warping
gips_project Landsat $ARGS --outdir landsat_project_nowarp --notld
gips_stats landsat_project_nowarp

# warp tiles
gips_tiles Landsat $ARGS --outdir landsat_warped_tiles --notld
gips_stats landsat_warped_tiles/*

# copy tiles
gips_tiles Landsat -t 012030 $DATES -v 4 --outdir landsat_tiles --notld
gips_stats landsat_tiles
