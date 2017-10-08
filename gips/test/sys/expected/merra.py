"""Known-good outcomes for tests, mostly stdout and created files."""

t_info = { 'stdout':  u"""\x1b[1mGIPS Data Repositories (v0.8.2)\x1b[0m
\x1b[1m
merra Products v1.0.0\x1b[0m
\x1b[1m
Standard Products
\x1b[0m   frland      Fraction of land (fraction)             
   patm        Surface atmospheric pressure (mb)       
   prcp        Daily total precipitation (mm day-1)    
   rhum        Relative humidity (%)                   
   shum        Relative humidity (kg kg-1)             
   srad        Incident solar radiation (W m-2)        
   tave        Ave daily air temperature data (K)      
   tmax        Max daily air temperature data          
   tmin        Min daily air temperature data          
   wind        Daily mean wind speed (m s-1)           
""" }


t_inventory = {
    'stdout': u"""\x1b[1mGIPS Data Inventory (v0.8.2)\x1b[0m
Retrieving inventory for site NHseacoast-0

\x1b[1mAsset Coverage for site NHseacoast-0\x1b[0m
\x1b[1m
Tile Coverage
\x1b[4m  Tile      % Coverage   % Tile Used\x1b[0m
  h01v01      100.0%        0.0%

\x1b[1m\x1b[4m    DATE       ASM       FLX       RAD       SLV     Product  \x1b[0m
\x1b[1m2015        
\x1b[0m    135                 100.0%     100.0%     100.0%   


1 files on 1 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mmerra: Modern Era Retrospective-analysis for Research and Applications\x1b[0m
"""
}


t_process = {
    'updated': {
        'merra/stage': None,
        'merra/tiles/h01v01/2015135': None,
    },
    'created': {
        'merra/tiles/h01v01/2015135/h01v01_2015135_merra_patm.tif': -832862627,
        'merra/tiles/h01v01/2015135/h01v01_2015135_merra_prcp.tif': -1814551991,
        'merra/tiles/h01v01/2015135/h01v01_2015135_merra_rhum.tif': 1203124684,
        'merra/tiles/h01v01/2015135/h01v01_2015135_merra_shum.tif': -97846975,
        'merra/tiles/h01v01/2015135/h01v01_2015135_merra_srad.tif': -1745572576,
        'merra/tiles/h01v01/2015135/h01v01_2015135_merra_tave.tif': 165650774,
        'merra/tiles/h01v01/2015135/h01v01_2015135_merra_tmax.tif': 463728356,
        'merra/tiles/h01v01/2015135/h01v01_2015135_merra_tmin.tif': -1469140067,
        'merra/tiles/h01v01/2015135/h01v01_2015135_merra_wind.tif': 1771906367,
    },
    'ignored': [
        'gips-inv-db.sqlite3',
    ],
    '_inv_stdout': """\x1b[1mGIPS Data Inventory (v0.8.2)\x1b[0m
Retrieving inventory for site NHseacoast-0

\x1b[1mAsset Coverage for site NHseacoast-0\x1b[0m
\x1b[1m
Tile Coverage
\x1b[4m  Tile      % Coverage   % Tile Used\x1b[0m
  h01v01      100.0%        0.0%

\x1b[1m\x1b[4m    DATE       ASM       FLX       RAD       SLV     Product  \x1b[0m
\x1b[1m2015        
\x1b[0m    135                 100.0%     100.0%     100.0%     \x1b[35mpatm\x1b[0m  \x1b[35mprcp\x1b[0m  \x1b[35mrhum\x1b[0m  \x1b[35mshum\x1b[0m  \x1b[35msrad\x1b[0m  \x1b[35mtave\x1b[0m  \x1b[35mtmax\x1b[0m  \x1b[35mtmin\x1b[0m  \x1b[35mwind\x1b[0m


1 files on 1 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mmerra: Modern Era Retrospective-analysis for Research and Applications\x1b[0m
""",
}

t_project = {
    'created': {
        '0': None,
        '0/2015135_merra_patm.tif': 147353283,
        '0/2015135_merra_prcp.tif': 1755012780,
        '0/2015135_merra_rhum.tif': -157842459,
        '0/2015135_merra_shum.tif': -1248597053,
        '0/2015135_merra_srad.tif': 335947606,
        '0/2015135_merra_tave.tif': 362688260,
        '0/2015135_merra_tmax.tif': -711449232,
        '0/2015135_merra_tmin.tif': 1704871784,
        '0/2015135_merra_wind.tif': 1730702417,
    }
}

t_project_no_warp = {
    'created': {
        '0': None,
        '0/2015135_merra_patm.tif': -1358396856,
        '0/2015135_merra_prcp.tif': -853804546,
        '0/2015135_merra_rhum.tif': 1810084110,
        '0/2015135_merra_shum.tif': -451670363,
        '0/2015135_merra_srad.tif': -872759650,
        '0/2015135_merra_tave.tif': 1065221007,
        '0/2015135_merra_tmax.tif': 1169791704,
        '0/2015135_merra_tmin.tif': 892785210,
        '0/2015135_merra_wind.tif': -15824969,
    }
}

# haven't ever used tiles
t_tiles = {
    'created': {
        'h01v01': None,
   },
}

t_tiles_copy = {
    'created': {
        'h01v01': None,
        'h01v01/h01v01_2015135_merra_patm.tif': -700661020,
        'h01v01/h01v01_2015135_merra_prcp.tif': -61013249,
        'h01v01/h01v01_2015135_merra_rhum.tif': 191746518,
        'h01v01/h01v01_2015135_merra_shum.tif': -1653055419,
        'h01v01/h01v01_2015135_merra_srad.tif': 777698408,
        'h01v01/h01v01_2015135_merra_tave.tif': 479810606,
        'h01v01/h01v01_2015135_merra_tmax.tif': 366977711,
        'h01v01/h01v01_2015135_merra_tmin.tif': 882403888,
        'h01v01/h01v01_2015135_merra_wind.tif': 301026612,
    }
}

t_stats = {
    'created': {
        'patm_stats.txt': 1367728250,
        'prcp_stats.txt': 1367728250,
        'rhum_stats.txt': 1367728250,
        'shum_stats.txt': 1367728250,
        'srad_stats.txt': 1367728250,
        'tave_stats.txt': 1367728250,
        'tmax_stats.txt': 1367728250,
        'tmin_stats.txt': 1367728250,
        'wind_stats.txt': 1367728250,
    }
}
