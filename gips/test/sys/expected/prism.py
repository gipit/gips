"""Known-good outcomes for tests, mostly stdout and created files."""

t_inventory = {
    'stdout': u"""\x1b[1mGIPS Data Inventory (v0.8.2)\x1b[0m
Retrieving inventory for site NHseacoast-0

\x1b[1mAsset Coverage for site NHseacoast-0\x1b[0m
\x1b[1m
Tile Coverage
\x1b[4m  Tile      % Coverage   % Tile Used\x1b[0m
   CONUS       99.2%        0.0%

\x1b[1m\x1b[4m    DATE       _ppt     _tmax     _tmin    Product  \x1b[0m
\x1b[1m1982        
\x1b[0m    335       99.2%     99.2%             
    336       99.2%     99.2%             
    337       99.2%     99.2%             


3 files on 3 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mprism: Daily Gridded Climate Data\x1b[0m
"""
}


t_process = {
    'updated': {
        'prism/tiles/CONUS/19821201': None,
        'prism/tiles/CONUS/19821202': None,
        'prism/tiles/CONUS/19821203': None
    },
    'created': {
        'prism/tiles/CONUS/19821201/CONUS_19821201_prism_ppt.tif': None,
        'prism/tiles/CONUS/19821201/CONUS_19821201_prism_tmax.tif': None,
        'prism/tiles/CONUS/19821201/PRISM_ppt_stable_4kmD2_19821201_bil.zip.index': 1277168171,
        'prism/tiles/CONUS/19821201/PRISM_tmax_stable_4kmD1_19821201_bil.zip.index': 267389174,
        'prism/tiles/CONUS/19821202/CONUS_19821202_prism_ppt.tif': None,
        'prism/tiles/CONUS/19821202/CONUS_19821202_prism_tmax.tif': None,
        'prism/tiles/CONUS/19821202/PRISM_ppt_stable_4kmD2_19821202_bil.zip.index': 1894168820,
        'prism/tiles/CONUS/19821202/PRISM_tmax_stable_4kmD1_19821202_bil.zip.index': -1880014768,
        'prism/tiles/CONUS/19821203/CONUS_19821203_prism_ppt.tif': None,
        'prism/tiles/CONUS/19821203/CONUS_19821203_prism_tmax.tif': None,
        'prism/tiles/CONUS/19821203/PRISM_ppt_stable_4kmD2_19821203_bil.zip.index': -764138626,
        'prism/tiles/CONUS/19821203/PRISM_tmax_stable_4kmD1_19821203_bil.zip.index': -267487182,
        'prism/tiles/CONUS/19821203/CONUS_19821203_prism_pptsum-3.tif': -1313341182
    },
}

