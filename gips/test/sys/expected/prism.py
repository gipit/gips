"""Known-good outcomes for tests, mostly stdout and created files."""

t_info = { 'stdout':  u"""\x1b[1mGIPS Data Repositories (v0.8.2)\x1b[0m
\x1b[1m
PRISM Products v0.1.1\x1b[0m
  Optional qualifiers listed below each product.
  Specify by appending '-option' to product (e.g., ref-toa)
\x1b[1m
Standard Products
\x1b[0m   ppt         Precipitate                             
   pptsum      Cumulative Precipitate                  
                 days: temporal window width (default: 3 days) .
   tmax        Daily Maximum Temperature               
   tmin        Daily Minimum Temperature               
""" }


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
\x1b[0m    335       99.2%     99.2%     99.2%   
    336       99.2%     99.2%     99.2%   
    337       99.2%     99.2%     99.2%   


3 files on 3 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mprism: Daily Gridded Climate Data\x1b[0m
""",
}


t_process = {
    'updated': {
        'prism/stage': None,
        'prism/tiles/CONUS/19821201': None,
        'prism/tiles/CONUS/19821202': None,
        'prism/tiles/CONUS/19821203': None
    },
    'created': {
        'prism/tiles/CONUS/19821201/CONUS_19821201_prism_ppt.tif': None,
        'prism/tiles/CONUS/19821201/CONUS_19821201_prism_tmax.tif': None,
        'prism/tiles/CONUS/19821201/CONUS_19821201_prism_tmin.tif': None,
        'prism/tiles/CONUS/19821202/CONUS_19821202_prism_ppt.tif': None,
        'prism/tiles/CONUS/19821202/CONUS_19821202_prism_tmax.tif': None,
        'prism/tiles/CONUS/19821202/CONUS_19821202_prism_tmin.tif': None,
        'prism/tiles/CONUS/19821203/CONUS_19821203_prism_ppt.tif': None,
        'prism/tiles/CONUS/19821203/CONUS_19821203_prism_tmax.tif': None,
        'prism/tiles/CONUS/19821203/CONUS_19821203_prism_tmin.tif': None,
        'prism/tiles/CONUS/19821203/CONUS_19821203_prism_pptsum-3.tif': 1770768528,
        'prism/tiles/CONUS/19821201/PRISM_tmax_stable_4kmD1_19821201_bil.zip.index': -74851230,
        'prism/tiles/CONUS/19821203/PRISM_ppt_stable_4kmD2_19821203_bil.zip.index': 1520480318,
        'prism/tiles/CONUS/19821203/PRISM_tmin_stable_4kmD1_19821203_bil.zip.index': -946867120,
        'prism/tiles/CONUS/19821202/PRISM_tmax_stable_4kmD1_19821202_bil.zip.index': 280913761,
        'prism/tiles/CONUS/19821202/PRISM_tmin_stable_4kmD1_19821202_bil.zip.index': 676175726,
        'prism/tiles/CONUS/19821203/PRISM_tmax_stable_4kmD1_19821203_bil.zip.index': -8647585,
        'prism/tiles/CONUS/19821202/PRISM_ppt_stable_4kmD2_19821202_bil.zip.index': -1251812608,
        'prism/tiles/CONUS/19821201/PRISM_tmin_stable_4kmD1_19821201_bil.zip.index': -1015405459,
        'prism/tiles/CONUS/19821201/PRISM_ppt_stable_4kmD2_19821201_bil.zip.index': 1582653443,
    },
    '_symlinks': {
        # TODO do the rest of the symlinks as time allows
        'prism/tiles/CONUS/19821201/CONUS_19821201_prism_ppt.tif':
            'prism/tiles/CONUS/19821201/PRISM_ppt_stable_4kmD2_19821201_bil.zip'
            '/PRISM_ppt_stable_4kmD2_19821201_bil.bil',
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
   CONUS       99.2%        0.0%

\x1b[1m\x1b[4m    DATE       _ppt     _tmax     _tmin    Product  \x1b[0m
\x1b[1m1982        
\x1b[0m    335       99.2%     99.2%     99.2%     \x1b[35mppt\x1b[0m  \x1b[35mtmax\x1b[0m  \x1b[35mtmin\x1b[0m
    336       99.2%     99.2%     99.2%     \x1b[35mppt\x1b[0m  \x1b[35mtmax\x1b[0m  \x1b[35mtmin\x1b[0m
    337       99.2%     99.2%     99.2%     \x1b[35mppt\x1b[0m  \x1b[35mpptsum-3\x1b[0m  \x1b[35mtmax\x1b[0m  \x1b[35mtmin\x1b[0m


3 files on 3 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mprism: Daily Gridded Climate Data\x1b[0m
""",
}

t_project = {
    'created': {
        '0': None,
        '0/1982335_prism_ppt.tif': -393322682,
        '0/1982335_prism_tmax.tif': 455141340,
        '0/1982335_prism_tmin.tif': 847249412,
        '0/1982336_prism_ppt.tif': 233388192,
        '0/1982336_prism_tmax.tif': 858443300,
        '0/1982336_prism_tmin.tif': -1078834087,
        '0/1982337_prism_ppt.tif': 983840913,
        '0/1982337_prism_pptsum.tif': 366592820,
        '0/1982337_prism_tmax.tif': 1878755022,
        '0/1982337_prism_tmin.tif': -100140734,
#        '0/1982335_prism_ppt.tif': -1525304978,
#        '0/1982335_prism_tmax.tif': 1259010415,
#        '0/1982335_prism_tmin.tif': -106507302,
#        '0/1982336_prism_ppt.tif': -92502584,
#        '0/1982336_prism_tmax.tif': 381286493,
#        '0/1982336_prism_tmin.tif': -801866142,
#        '0/1982337_prism_ppt.tif': 670927499,
#        '0/1982337_prism_pptsum.tif': 1195812568,
#        '0/1982337_prism_tmax.tif': 1464461372,
#        '0/1982337_prism_tmin.tif': -2082314142,
    }
}

t_project_no_warp = {
    'created': {'0': None,
        '0/1982335_prism_ppt.tif': -1700182526,
        '0/1982337_prism_tmin.tif': 2101752746,
        '0/1982336_prism_ppt.tif': -1700182526,
        '0/1982337_prism_ppt.tif': -1352376509,
        '0/1982335_prism_tmax.tif': 1306198481,
        '0/1982336_prism_tmax.tif': 1773055017,
        '0/1982337_prism_tmax.tif': -1029779846,
        '0/1982335_prism_tmin.tif': 976342700,
        '0/1982336_prism_tmin.tif': -165045365,
        '0/1982337_prism_pptsum.tif': -1460879448,
#        '0/1982335_prism_ppt.tif': -427336450,
#        '0/1982335_prism_tmax.tif': 571279318,
#        '0/1982335_prism_tmin.tif': 1498899978,
#        '0/1982336_prism_ppt.tif': -427336450,
#        '0/1982336_prism_tmax.tif': -1198323967,
#        '0/1982336_prism_tmin.tif': 88143743,
#        '0/1982337_prism_ppt.tif': -1329047557,
#        '0/1982337_prism_pptsum.tif': -1792182993,
#        '0/1982337_prism_tmax.tif': -622113590,
#        '0/1982337_prism_tmin.tif': -1048138502,
    }
}

# haven't ever used tiles
# t_tiles = {
#     'created': {},
# }

t_tiles_copy = {
    'created': {
        'CONUS': None,
        'CONUS/CONUS_1982335_prism_ppt.tif': -55072421,
        'CONUS/CONUS_1982335_prism_tmax.tif': -234123443,
        'CONUS/CONUS_1982335_prism_tmin.tif': -1003055785,
        'CONUS/CONUS_1982336_prism_ppt.tif': -1365707606,
        'CONUS/CONUS_1982336_prism_tmax.tif': 2008098128,
        'CONUS/CONUS_1982336_prism_tmin.tif': -1569456782,
        'CONUS/CONUS_1982337_prism_ppt.tif': 570845567,
        'CONUS/CONUS_1982337_prism_tmax.tif': 1289336243,
        'CONUS/CONUS_1982337_prism_tmin.tif': -818489243,
        'CONUS/CONUS_1982337_prism_pptsum.tif': 272867928,
    }
}

t_stats = { 
    'created': {
        'ppt_stats.txt': -2108420090,
        'pptsum_stats.txt': -1221703895,
        'tmax_stats.txt': -628985934,
        'tmin_stats.txt': 691861124
    }
}

