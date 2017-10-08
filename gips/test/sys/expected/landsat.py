"""Known-good outcomes for tests, mostly stdout and created files."""

# trailing whitespace and other junk characters are in current output
t_info = { 'stdout':  u"""\x1b[1mGIPS Data Repositories (v0.8.2)\x1b[0m
\x1b[1m
Landsat Products v1.0.0\x1b[0m
  Optional qualifiers listed below each product.
  Specify by appending '-option' to product (e.g., ref-toa)
\x1b[1m
Index Products
\x1b[0m   bi          Brightness Index                        
                 toa: use top of the atmosphere reflectance
   evi         Enhanced Vegetation Index               
                 toa: use top of the atmosphere reflectance
   lswi        Land Surface Water Index                
                 toa: use top of the atmosphere reflectance
   msavi2      Modified Soil-Adjusted Vegetation Index (revised)
                 toa: use top of the atmosphere reflectance
   ndsi        Normalized Difference Snow Index        
                 toa: use top of the atmosphere reflectance
   ndvi        Normalized Difference Vegetation Index  
                 toa: use top of the atmosphere reflectance
   ndwi        Normalized Difference Water Index       
                 toa: use top of the atmosphere reflectance
   satvi       Soil-Adjusted Total Vegetation Index    
                 toa: use top of the atmosphere reflectance
   vari        Visible Atmospherically Resistant Index 
                 toa: use top of the atmosphere reflectance
\x1b[1m
LC8SR Products
\x1b[0m   ndvi8sr     Normalized Difference Vegetation from LC8SR
\x1b[1m
ACOLITE Products
\x1b[0m   acoflags    0 = water 1 = no data 2 = land          
   fai         Floating Algae Index                    
   oc2chl      Blue-Green Ratio Chlorophyll Algorithm using bands 483 & 561
   oc3chl      Blue-Green Ratio Chlorophyll Algorithm using bands 443, 483, & 561
   rhow        Water-Leaving Radiance-Reflectance      
   spm655      Suspended Sediment Concentration 655nm  
   turbidity   Blended Turbidity                       
\x1b[1m
Tillage Products
\x1b[0m   crc         Crop Residue Cover                      
                 toa: use top of the atmosphere reflectance
   isti        Inverse Standard Tillage Index          
                 toa: use top of the atmosphere reflectance
   ndti        Normalized Difference Tillage Index     
                 toa: use top of the atmosphere reflectance
   sti         Standard Tillage Index                  
                 toa: use top of the atmosphere reflectance
\x1b[1m
Standard Products
\x1b[0m   acca        Automated Cloud Cover Assessment        
                 X: erosion kernel diameter in pixels (default: 5)
                 Y: dilation kernel diameter in pixels (default: 10)
                 Z: cloud height in meters (default: 4000)
   bqa         The quality band extracted into separate layers.
   bqashadow   LC8 QA + Shadow Smear                   
                 X: erosion kernel diameter in pixels (default: 5)
                 Y: dilation kernel diameter in pixels (default: 10)
                 Z: cloud height in meters (default: 4000)
   cloudmask   Cloud (and shadow) mask product based on cloud bits of the quality band
   dn          Raw digital numbers                     
   fmask       Fmask cloud cover                       
   landmask    Land mask from LC8SR                    
   rad         Surface-leaving radiance                
                 toa: use top of the atmosphere reflectance
   ref         Surface reflectance                     
                 toa: use top of the atmosphere reflectance
   tcap        Tassled cap transformation              
   temp        Brightness (apparent) temperature       
   volref      Volumetric water reflectance - valid for water only
                 toa: use top of the atmosphere reflectance
   wtemp       Water temperature (atmospherically correct) - valid for water only
"""}

t_inventory = { 'stdout': u"""\x1b[1mGIPS Data Inventory (v0.8.2)\x1b[0m
Retrieving inventory for site NHseacoast-0
fname
LC08_L1GT_012030_20151218_20170224_01_T2.tar.gz
C1 asset

\x1b[1mAsset Coverage for site NHseacoast-0\x1b[0m
\x1b[1m
Tile Coverage
\x1b[4m  Tile      % Coverage   % Tile Used\x1b[0m
  012030      100.0%        6.7%
  013030        2.4%        0.2%

\x1b[1m\x1b[4m    DATE        C1        DN        SR     Product  \x1b[0m
\x1b[1m2015        
\x1b[0m    352       100.0%                       


1 files on 1 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mLC8: Landsat 8\x1b[0m
\x1b[31mLC8SR: Landsat 8 Surface Reflectance\x1b[0m
\x1b[32mLE7: Landsat 7\x1b[0m
\x1b[34mLT5: Landsat 5\x1b[0m
"""}

t_process = {
    'compare_stderr': False,
    'updated': {
        'landsat/stage': None,
        'landsat/tiles/012030/2015352': None
    },
    'created': {
        'landsat/tiles/012030/2015352/012030_2015352_LC8_acca.tif': -531492048,
        'landsat/tiles/012030/2015352/012030_2015352_LC8_bqashadow.tif': -1819149482,
        'landsat/tiles/012030/2015352/012030_2015352_LC8_ndvi-toa.tif': 329107382,
        'landsat/tiles/012030/2015352/012030_2015352_LC8_rad-toa.tif': -1222249885,
        'landsat/tiles/012030/2015352/012030_2015352_LC8_ref-toa.tif': -871936054,
        'landsat/tiles/012030/2015352/LC08_L1GT_012030_20151218_20170224_01_T2.tar.gz.index': -394988487,
        'landsat/tiles/012030/2015352/LC08_L1GT_012030_20151218_20170224_01_T2_MTL.txt': -1453474890,
    },
    'ignored': [
        'gips-inv-db.sqlite3',
    ]
}

t_process_acolite = {
    'created': {
        'landsat/tiles/012030/2017213/012030_2017213_LC8_acoflags.tif': -514981863,
        'landsat/tiles/012030/2017213/012030_2017213_LC8_fai.tif': -1820016301,
        'landsat/tiles/012030/2017213/012030_2017213_LC8_oc2chl.tif': -794776169,
        'landsat/tiles/012030/2017213/012030_2017213_LC8_oc3chl.tif': 337875725,
        'landsat/tiles/012030/2017213/012030_2017213_LC8_rhow.tif': 1044772611,
        'landsat/tiles/012030/2017213/012030_2017213_LC8_spm655.tif': 1437263134,
        'landsat/tiles/012030/2017213/012030_2017213_LC8_turbidity.tif': -1347918033,
    },
    'updated': {
        'landsat/stage': None,
        'landsat/tiles/012030/2017213': None
    },
    'ignored': [
        'gips-inv-db.sqlite3',
    ],
}

t_project = {
    'compare_stderr': False,
    'created': {
        '0': None,
        '0/2015352_LC8_acca.tif': 402348046,
        '0/2015352_LC8_bqashadow.tif': 923940030,
        '0/2015352_LC8_ndvi-toa.tif': 728893178,
        '0/2015352_LC8_rad-toa.tif': -1053542955,
        '0/2015352_LC8_ref-toa.tif': -1149010214,
    }
}

t_project_no_warp = {
    'compare_stderr': False,
    'created': {
        '0': None,
        '0/2015352_LC8_acca.tif': -126711306,
        '0/2015352_LC8_bqashadow.tif': 1681911857,
        '0/2015352_LC8_ndvi-toa.tif': 1662486138,
        '0/2015352_LC8_rad-toa.tif': -196115636,
        '0/2015352_LC8_ref-toa.tif': -1147999741,
    }
}

# TODO this bug rearing its ugly head again?
# See https://github.com/Applied-GeoSolutions/gips/issues/54
t_tiles = { 'created': {'012030': None}}

t_tiles_copy = {
    'compare_stderr': False,
    'created': {
        '012030': None,
        '012030/012030_2015352_LC8_acca.tif': 176561467,
        '012030/012030_2015352_LC8_bqashadow.tif': 912021217,
        '012030/012030_2015352_LC8_ndvi-toa.tif': -1333295744,
        '012030/012030_2015352_LC8_rad-toa.tif': 1609412102,
        '012030/012030_2015352_LC8_ref-toa.tif': -1797834447,
    }
}

t_stats = { 'created': {
    'acca_stats.txt': -174967201,
    'bqashadow_stats.txt': 1868908586,
    'ndvi-toa_stats.txt': -1084861813,
    'rad-toa_stats.txt': -545320378,
    'ref-toa_stats.txt': -1132928652,
}}
