"""Known-good outcomes for tests, mostly stdout and created files."""

# trailing whitespace and other junk characters are in current output
t_info = { 'stdout':  """\x1b[1mGIPS Data Repositories (v0.8.2)\x1b[0m
\x1b[1m
Landsat Products v0.9.0\x1b[0m
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
\x1b[1m
LC8SR Products
\x1b[0m   ndvi8sr     Normalized Difference Vegetation from LC8SR
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
   bqa         LC8 band quality                        
   bqashadow   LC8 QA + Shadow Smear                   
                 X: erosion kernel diameter in pixels (default: 5)
                 Y: dilation kernel diameter in pixels (default: 10)
                 Z: cloud height in meters (default: 4000)
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

t_inventory = { 'stdout': """\x1b[1mGIPS Data Inventory (v0.8.2)\x1b[0m
Retrieving inventory for site NHseacoast-0
fname
LC80120302015352LGN00.tar.gz
DN asset

\x1b[1mAsset Coverage for site NHseacoast-0\x1b[0m
\x1b[1m
Tile Coverage
\x1b[4m  Tile      % Coverage   % Tile Used\x1b[0m
  012030      100.0%        6.7%
  013030        2.4%        0.2%

\x1b[1m\x1b[4m    DATE        DN        SR     Product  \x1b[0m
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
    'updated': {'landsat/tiles/012030/2015352': None},
    'created': {
        'landsat/tiles/012030/2015352/012030_2015352_LC8_acca.tif': -1532119925,
        'landsat/tiles/012030/2015352/012030_2015352_LC8_bqashadow.tif': -1566662956,
        'landsat/tiles/012030/2015352/012030_2015352_LC8_ndvi-toa.tif': -204896503,
        'landsat/tiles/012030/2015352/012030_2015352_LC8_rad-toa.tif': 1658217796,
        'landsat/tiles/012030/2015352/012030_2015352_LC8_ref-toa.tif': -1766137404,
        'landsat/tiles/012030/2015352/LC80120302015352LGN00.tar.gz.index': 1896263933,
        'landsat/tiles/012030/2015352/LC80120302015352LGN00_MTL.txt': 2084350553,
    }
}

t_project = {
    'compare_stderr': False,
    'created': {
        '0': None,
        '0/2015352_LC8_acca.tif': -1824905460,
        '0/2015352_LC8_bqashadow.tif': 1603304372,
        '0/2015352_LC8_ndvi-toa.tif': 844246796,
        '0/2015352_LC8_rad-toa.tif': -317896577,
        '0/2015352_LC8_ref-toa.tif': -1496156246,
    }
}

t_project_no_warp = {
    'compare_stderr': False,
    'created': {
        '0': None,
        '0/2015352_LC8_acca.tif': 1288527028,
        '0/2015352_LC8_bqashadow.tif': -2134364541,
        '0/2015352_LC8_ndvi-toa.tif': 1466287813,
        '0/2015352_LC8_rad-toa.tif': 1618019596,
        '0/2015352_LC8_ref-toa.tif': -614392204,
    }
}

# TODO this bug rearing its ugly head again?
# See https://github.com/Applied-GeoSolutions/gips/issues/54
t_tiles = { 'created': {'012030': None}}

t_tiles_copy = {
    'compare_stderr': False,
    'created': {
        '012030': None,
	'012030/012030_2015352_LC8_acca.tif': 1593997289,
	'012030/012030_2015352_LC8_bqashadow.tif': -116702962,
	'012030/012030_2015352_LC8_ndvi-toa.tif': -1923562909,
	'012030/012030_2015352_LC8_rad-toa.tif': 525640715,
	'012030/012030_2015352_LC8_ref-toa.tif': -1287648712,
    }
}

t_stats = { 'created': {
    'acca_stats.txt': -789655715,
    'bqashadow_stats.txt': 1501756012,
    'ndvi-toa_stats.txt': -77721729,
    'rad-toa_stats.txt': 1664250177,
    'ref-toa_stats.txt': 2007199405,
}}
