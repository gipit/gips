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

_en_base = ('sentinel2/tiles/19TCH/2017010'
            '/S2A_MSIL1C_20170110T153611_N0204_R111_T19TCH_20170110T154029.SAFE')

_en_prefix = _en_base + '/GRANULE/L1C_T19TCH_A008116_20170110T154029/IMG_DATA'

_extraction_noise = [
    _en_base,
    _en_base + '/GRANULE',
    _en_base + '/GRANULE/L1C_T19TCH_A008116_20170110T154029',
    _en_prefix,
    _en_prefix + '/T19TCH_20170110T153611_B01.jp2',
    _en_prefix + '/T19TCH_20170110T153611_B02.jp2',
    _en_prefix + '/T19TCH_20170110T153611_B02.jp2.aux.xml',
    _en_prefix + '/T19TCH_20170110T153611_B03.jp2',
    _en_prefix + '/T19TCH_20170110T153611_B03.jp2.aux.xml',
    _en_prefix + '/T19TCH_20170110T153611_B04.jp2',
    _en_prefix + '/T19TCH_20170110T153611_B04.jp2.aux.xml',
    _en_prefix + '/T19TCH_20170110T153611_B05.jp2',
    _en_prefix + '/T19TCH_20170110T153611_B06.jp2',
    _en_prefix + '/T19TCH_20170110T153611_B07.jp2',
    _en_prefix + '/T19TCH_20170110T153611_B08.jp2',
    _en_prefix + '/T19TCH_20170110T153611_B08.jp2.aux.xml',
    _en_prefix + '/T19TCH_20170110T153611_B09.jp2',
    _en_prefix + '/T19TCH_20170110T153611_B10.jp2',
    _en_prefix + '/T19TCH_20170110T153611_B11.jp2',
    _en_prefix + '/T19TCH_20170110T153611_B11.jp2.aux.xml',
    _en_prefix + '/T19TCH_20170110T153611_B12.jp2',
    _en_prefix + '/T19TCH_20170110T153611_B12.jp2.aux.xml',
    _en_prefix + '/T19TCH_20170110T153611_B8A.jp2',
]

t_process = {
    'compare_stderr': False,
    'updated': {
        'sentinel2/stage': None,
        'sentinel2/tiles/19TCH/2017010': None,
    },
    'created': {
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
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_ref-toa.tif': -375067000,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_satvi-toa.tif': 1636524188,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_sti-toa.tif': 530354401,
        'sentinel2/tiles/19TCH/2017010/19TCH_2017010_S2A_vari-toa.tif': 1830182386,
    },
    'ignored': _extraction_noise + [
        'gips-inv-db.sqlite3',
        'sentinel2/tiles/19TCH/2017010/S2A_MSIL1C_20170110T153611_N0204_R111_T19TCH_20170110T154029.zip',
        'sentinel2/tiles/19TCH/2017010/S2A_MSIL1C_20170110T153611_N0204_R111_T19TCH_20170110T154029.zip.index',
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
   brgt        VIS and NIR reflectance, weighted by solar energy distribution.
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
\x1b[0m   rad         Surface-leaving radiance                
   ref         Surface reflectance                     
"""
}


t_project = {
    'compare_stderr': False,
    'created': {
        '0': None,
        '0/2017010_S2A_bi-toa.tif': -29689498,
        '0/2017010_S2A_brgt-toa.tif': -772350197,
        '0/2017010_S2A_crc-toa.tif': -1762466646,
        '0/2017010_S2A_crcm-toa.tif': -1642853419,
        '0/2017010_S2A_evi-toa.tif': 898489308,
        '0/2017010_S2A_isti-toa.tif': -308268452,
        '0/2017010_S2A_lswi-toa.tif': 970713245,
        '0/2017010_S2A_msavi2-toa.tif': 1140376127,
        '0/2017010_S2A_ndsi-toa.tif': 692688448,
        '0/2017010_S2A_ndti-toa.tif': -1339243026,
        '0/2017010_S2A_ndvi-toa.tif': 290485324,
        '0/2017010_S2A_ref-toa.tif': 875075934,
        '0/2017010_S2A_satvi-toa.tif': -2092483419,
        '0/2017010_S2A_sti-toa.tif': -2123877142,
        '0/2017010_S2A_vari-toa.tif': -2127272542,
    },
    'ignored': _extraction_noise,
}


t_stats = {
    'created': {
        'bi-toa_stats.txt': 182802308,
        'brgt-toa_stats.txt': 182802308,
        'crc-toa_stats.txt': -1227739799,
        'crcm-toa_stats.txt': 1974395793,
        'evi-toa_stats.txt': -1853135550,
        'isti-toa_stats.txt': -95638900,
        'lswi-toa_stats.txt': 1301066971,
        'msavi2-toa_stats.txt': 1963970515,
        'ndsi-toa_stats.txt': 749844046,
        'ndti-toa_stats.txt': -1079013241,
        'ndvi-toa_stats.txt': -1186889206,
        'ref-toa_stats.txt': 711244685,
        'satvi-toa_stats.txt': 154980132,
        'sti-toa_stats.txt': 182802308,
        'vari-toa_stats.txt': -462847554,
    },
    'ignored': _extraction_noise,
}

t_process_acolite = {
    'compare_stderr': False,
    'updated': {
        'sentinel2/stage': None,
        'sentinel2/tiles/19TCH/2017010': None,
    },
    'created': {
        # TODO
    },
    'ignored': _extraction_noise + [
        'gips-inv-db.sqlite3',
        'sentinel2/tiles/19TCH/2017010/S2A_MSIL1C_20170110T153611_N0204_R111_T19TCH_20170110T154029.zip',
        'sentinel2/tiles/19TCH/2017010/S2A_MSIL1C_20170110T153611_N0204_R111_T19TCH_20170110T154029.zip.index',
    ],
}
