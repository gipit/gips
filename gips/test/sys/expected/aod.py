"""Known-good outcomes for tests, mostly stdout and created files."""

t_info = { 'stdout':  u"""\x1b[1mGIPS Data Repositories (v0.8.2)\x1b[0m
\x1b[1m
AOD Products v0.9.1\x1b[0m
\x1b[1m
Standard Products
\x1b[0m   aod         Aerosol Optical Depth                   
   lta         Average AOD                             
   ltad        Average daily AOD                       
""" }


t_inventory = {
    'stdout': u"""\x1b[1mGIPS Data Inventory (v0.8.2)\x1b[0m
Retrieving inventory for site NE-0

\x1b[1mAsset Coverage for site NE-0\x1b[0m
\x1b[1m
Tile Coverage
\x1b[4m  Tile      % Coverage   % Tile Used\x1b[0m
  h01v01      100.0%      100.0%

\x1b[1m\x1b[4m    DATE      MOD08    Product  \x1b[0m
\x1b[1m2017        
\x1b[0m    004       100.0%     \x1b[35maod\x1b[0m
    005       100.0%     \x1b[35maod\x1b[0m
    006       100.0%     \x1b[35maod\x1b[0m


3 files on 3 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mMOD: MODIS Terra\x1b[0m
""",
}


t_process = {
    'updated': {
    },
    'created': {
    },
    'ignored': [
        'gips-inv-db.sqlite3',
    ],
    'compare_stdout': False,
    '_inv_stdout': """\x1b[1mGIPS Data Inventory (v0.8.2)\x1b[0m
Retrieving inventory for site NE-0

\x1b[1mAsset Coverage for site NE-0\x1b[0m
\x1b[1m
Tile Coverage
\x1b[4m  Tile      % Coverage   % Tile Used\x1b[0m
  h01v01      100.0%      100.0%

\x1b[1m\x1b[4m    DATE      MOD08    Product  \x1b[0m
\x1b[1m2017        
\x1b[0m    004       100.0%     \x1b[35maod\x1b[0m
    005       100.0%     \x1b[35maod\x1b[0m
    006       100.0%     \x1b[35maod\x1b[0m


3 files on 3 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mMOD: MODIS Terra\x1b[0m
""",
}

t_project = {
    'compare_stdout': False,
    'compare_stderr': False,
    'created': {
        '0': None,
        '0/2017004_MOD_aod.tif': -766418958,
        '0/2017005_MOD_aod.tif': 1283286056,
        '0/2017006_MOD_aod.tif': -2004417139,
    }
}

t_project_no_warp = {
    'compare_stdout': False,
    'created': {
        '0': None,
        '0/2017004_MOD_aod.tif': -124661770,
        '0/2017005_MOD_aod.tif': -1798929705,
        '0/2017006_MOD_aod.tif': 1354760932,
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
            'h01v01/h01v01_2017004_MOD_aod.tif': 190937103,
            'h01v01/h01v01_2017005_MOD_aod.tif': -1587512880,
            'h01v01/h01v01_2017006_MOD_aod.tif': -2112843720,
    }
}

t_stats = {
    'stdout': u"""\x1b[1mGIPS Image Statistics (v0.1.0)\x1b[0m
Stats for Project directory: /tmp/test-output
Calculating statistics for 2017-01-04
Calculating statistics for 2017-01-05
Calculating statistics for 2017-01-06
""",
    'created': {
        'aod_stats.txt': -769515142
    }
}
