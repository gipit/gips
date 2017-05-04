"""Known-good outcomes for tests, mostly stdout and created files."""

t_inventory = { 'stdout': u"""\x1b[1mGIPS Data Inventory (v0.8.2)\x1b[0m
Retrieving inventory for site durham-0

\x1b[1mAsset Coverage for site durham-0\x1b[0m
\x1b[1m
Tile Coverage
\x1b[4m  Tile      % Coverage   % Tile Used\x1b[0m
   19TCH      100.0%        0.1%

\x1b[1m\x1b[4m    DATE       L1C     Product  \x1b[0m
\x1b[1m2017        
\x1b[0m    010       100.0%   


1 files on 1 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mS2A: Sentinel-2, Satellite A\x1b[0m
\x1b[31mS2B: Sentinel-2, Satellite B\x1b[0m
"""
}

t_process = {
    'compare_stderr': False,
    'updated': {
        'sentinel2/tiles/19TCH/2017010': None,
    },
    'created': {
        # jan 10
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_bi-toa.tif': 541673910,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_brgt-toa.tif': 38425695,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_crc-toa.tif': -692365694,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_crcm-toa.tif': -838929383,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_evi-toa.tif': -850280958,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_isti-toa.tif': -1580129686,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_lswi-toa.tif': 1676471139,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_msavi2-toa.tif': -536630164,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_ndsi-toa.tif': -2118476257,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_ndti-toa.tif': -716676364,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_ndvi-toa.tif': 1460522898,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_ref-toa.tif': -711062552,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_satvi-toa.tif': 1636524188,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_sti-toa.tif': 530354401,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_vari-toa.tif': 1830182386,
        # jan 30
        #'sentinel2/tiles/19TCH/2017030/19TCH_2017030_S2A_bi-toa.tif': -1803674796,
        #'sentinel2/tiles/19TCH/2017030/19TCH_2017030_S2A_brgt-toa.tif': 38425695,
        #'sentinel2/tiles/19TCH/2017030/19TCH_2017030_S2A_crc-toa.tif': 1709338205,
        #'sentinel2/tiles/19TCH/2017030/19TCH_2017030_S2A_crcm-toa.tif': -409438165,
        #'sentinel2/tiles/19TCH/2017030/19TCH_2017030_S2A_evi-toa.tif': 323913491,
        #'sentinel2/tiles/19TCH/2017030/19TCH_2017030_S2A_isti-toa.tif': -1749541518,
        #'sentinel2/tiles/19TCH/2017030/19TCH_2017030_S2A_lswi-toa.tif': 2114258079,
        #'sentinel2/tiles/19TCH/2017030/19TCH_2017030_S2A_msavi2-toa.tif': 2091111873,
        #'sentinel2/tiles/19TCH/2017030/19TCH_2017030_S2A_ndsi-toa.tif': 1825007655,
        #'sentinel2/tiles/19TCH/2017030/19TCH_2017030_S2A_ndti-toa.tif': 42954357,
        #'sentinel2/tiles/19TCH/2017030/19TCH_2017030_S2A_ndvi-toa.tif': 1854221224,
        #'sentinel2/tiles/19TCH/2017030/19TCH_2017030_S2A_ref-toa.tif': -1762824208,
        #'sentinel2/tiles/19TCH/2017030/19TCH_2017030_S2A_satvi-toa.tif': 102386778,
        #'sentinel2/tiles/19TCH/2017030/19TCH_2017030_S2A_sti-toa.tif': 1625215040,
        #'sentinel2/tiles/19TCH/2017030/19TCH_2017030_S2A_vari-toa.tif': -1662958906,
    },
    'ignored': [
        'gips-inv-db.sqlite3',
        'sentinel2/tiles/19TCH/2017010/S2A_MSIL1C_20170110T153611_N0204_R111_T19TCH_20170110T154029.zip',
        'sentinel2/tiles/19TCH/2017030',
        'sentinel2/tiles/19TCH/2017030/S2A_MSIL1C_20170130T153551_N0204_R111_T19TCH_20170130T153712.zip',
    ],
    '_inv_stdout': u"""\x1b[1mGIPS Data Inventory (v0.8.2)\x1b[0m
Retrieving inventory for site durham-0

\x1b[1mAsset Coverage for site durham-0\x1b[0m
\x1b[1m
Tile Coverage
\x1b[4m  Tile      % Coverage   % Tile Used\x1b[0m
   19TCH      100.0%        0.1%

\x1b[1m\x1b[4m    DATE       L1C     Product  \x1b[0m
\x1b[1m2017        
\x1b[0m    010       100.0%     \x1b[35mbi-toa\x1b[0m  \x1b[35mbrgt-toa\x1b[0m  \x1b[35mcrc-toa\x1b[0m  \x1b[35mcrcm-toa\x1b[0m  \x1b[35mevi-toa\x1b[0m  \x1b[35misti-toa\x1b[0m  \x1b[35mlswi-toa\x1b[0m  \x1b[35mmsavi2-toa\x1b[0m  \x1b[35mndsi-toa\x1b[0m  \x1b[35mndti-toa\x1b[0m  \x1b[35mndvi-toa\x1b[0m  \x1b[35mref-toa\x1b[0m  \x1b[35msatvi-toa\x1b[0m  \x1b[35msti-toa\x1b[0m  \x1b[35mvari-toa\x1b[0m


1 files on 1 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mS2A: Sentinel-2, Satellite A\x1b[0m
\x1b[31mS2B: Sentinel-2, Satellite B\x1b[0m
"""
}


t_info = { 'stdout': u"""\x1b[1mGIPS Data Repositories (v0.8.2)\x1b[0m
\x1b[1m
Sentinel2 Products v0.1.0\x1b[0m
\x1b[1m
Index Products
\x1b[0m   bi          Brightness Index                        
   brgt        Brightness index:  Visible to near infrared reflectance weighted by approximate energy distribution of the solar spectrum. A proxy for broadband albedo.
   crc         Crop Residue Cover (uses BLUE)          
   crcm        Crop Residue Cover, Modified (uses GREEN)
   evi         Enhanced Vegetation Index               
   isti        Inverse Standard Tillage Index          
   lswi        Land Surface Water Index                
   msavi2      Modified Soil-adjusted Vegetation Index 
   ndsi        Normalized Difference Snow Index        
   ndti        Normalized Difference Tillage Index     
   ndvi        Normalized Difference Vegetation Index  
   satvi       Soil-Adjusted Total Vegetation Index    
   sti         Standard Tillage Index                  
   vari        Visible Atmospherically Resistant Index 
\x1b[1m
Standard Products
\x1b[0m   ref         Surface reflectance                     
"""
}


t_project = {
    'created': {
        '0': None,
        '0/2017010_S2A_ref-toa.tif': 56081467,
        '0/2017010_S2A_bi-toa.tif': -531404086,
        '0/2017010_S2A_brgt-toa.tif': 1949224049,
        '0/2017010_S2A_crc-toa.tif': -1921980023,
        '0/2017010_S2A_crcm-toa.tif': -620088516,
        '0/2017010_S2A_evi-toa.tif': 1327644820,
        '0/2017010_S2A_isti-toa.tif': -621591635,
        '0/2017010_S2A_lswi-toa.tif': 1687749861,
        '0/2017010_S2A_msavi2-toa.tif': -1672549497,
        '0/2017010_S2A_ndsi-toa.tif': 1678838722,
        '0/2017010_S2A_ndti-toa.tif': -1245896923,
        '0/2017010_S2A_ndvi-toa.tif': 754895178,
        '0/2017010_S2A_satvi-toa.tif': -339325105,
        '0/2017010_S2A_sti-toa.tif': 528679666,
        '0/2017010_S2A_vari-toa.tif': -915656392,
    },
}


t_stats = { 'created': {
    'bi-toa_stats.txt': 1420978991,
    'brgt-toa_stats.txt': 1420978991,
    'crc-toa_stats.txt': 1420978991,
    'crcm-toa_stats.txt': 1420978991,
    'evi-toa_stats.txt': 1420978991,
    'isti-toa_stats.txt': 1420978991,
    'lswi-toa_stats.txt': 1420978991,
    'msavi2-toa_stats.txt': 1420978991,
    'ndsi-toa_stats.txt': 1420978991,
    'ndti-toa_stats.txt': 1420978991,
    'ndvi-toa_stats.txt': 1420978991,
    'ref-toa_stats.txt': 698328601,
    'satvi-toa_stats.txt': 1420978991,
    'sti-toa_stats.txt': 1420978991,
    'vari-toa_stats.txt': 1420978991,
}}
