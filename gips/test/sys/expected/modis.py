"""Known-good outcomes for tests, mostly stdout and created files."""

t_inventory = { 'stdout': u"""\x1b[1mGIPS Data Inventory (v0.8.2)\x1b[0m
Retrieving inventory for site NHseacoast-0

\x1b[1mAsset Coverage for site NHseacoast-0\x1b[0m
\x1b[1m
Tile Coverage
\x1b[4m  Tile      % Coverage   % Tile Used\x1b[0m
  h12v04      100.0%        0.2%

\x1b[1m\x1b[4m    DATE     MCD12Q1   MCD43A2   MCD43A4   MOD09Q1   MOD10A1   MOD10A2   MOD11A1   MOD11A2   MYD10A1   MYD10A2   MYD11A1   MYD11A2   Product  \x1b[0m
\x1b[1m2012        
\x1b[0m    336                 100.0%     100.0%               100.0%               100.0%               100.0%               100.0%             
    337                 100.0%     100.0%     100.0%     100.0%               100.0%     100.0%     100.0%               100.0%             
    338                 100.0%     100.0%               100.0%               100.0%               100.0%               100.0%             


3 files on 3 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mMCD: Aqua/Terra Combined\x1b[0m
\x1b[31mMOD: Terra\x1b[0m
\x1b[32mMOD-MYD: Aqua/Terra together\x1b[0m
\x1b[34mMYD: Aqua\x1b[0m
"""
}


t_process = {
    'updated': {'modis/tiles/h12v04/2012336': None,
             'modis/tiles/h12v04/2012337': None,
             'modis/tiles/h12v04/2012338': None},
    'created': {
        # TODO Are these broken or what?  Each None is a broken symlink that seems to rely on weird
        # gdal behavior (it uses the broken target of the symlink as a sort of data identifier?):
        # See https://github.com/Applied-GeoSolutions/gips/issues/54
        'modis/tiles/h12v04/2012336/h12v04_2012336_MCD_quality.tif': None,
        'modis/tiles/h12v04/2012337/h12v04_2012337_MCD_quality.tif': None,
        'modis/tiles/h12v04/2012337/h12v04_2012337_MOD_temp8td.tif': None,
        'modis/tiles/h12v04/2012337/h12v04_2012337_MOD_temp8tn.tif': None,
        'modis/tiles/h12v04/2012338/h12v04_2012338_MCD_quality.tif': None,
        'modis/tiles/h12v04/2012336/h12v04_2012336_MCD_fsnow.tif': -843500181,
        'modis/tiles/h12v04/2012336/h12v04_2012336_MCD_indices.tif': -684782438,
        'modis/tiles/h12v04/2012336/h12v04_2012336_MCD_snow.tif': 388495321,
        'modis/tiles/h12v04/2012336/h12v04_2012336_MOD-MYD_obstime.tif': -661225282,
        'modis/tiles/h12v04/2012336/h12v04_2012336_MOD-MYD_temp.tif': 1278807725,
        'modis/tiles/h12v04/2012336/h12v04_2012336_MOD_clouds.tif': 161070470,
        'modis/tiles/h12v04/2012337/h12v04_2012337_MCD_fsnow.tif': 297883486,
        'modis/tiles/h12v04/2012337/h12v04_2012337_MCD_indices.tif': -1806386481,
        'modis/tiles/h12v04/2012337/h12v04_2012337_MCD_snow.tif': -748640537,
        'modis/tiles/h12v04/2012337/h12v04_2012337_MOD-MYD_obstime.tif': 348598487,
        'modis/tiles/h12v04/2012337/h12v04_2012337_MOD-MYD_temp.tif': 1853031589,
        'modis/tiles/h12v04/2012337/h12v04_2012337_MOD_clouds.tif': -832284681,
        'modis/tiles/h12v04/2012337/h12v04_2012337_MOD_ndvi8.tif': -395723296,
        'modis/tiles/h12v04/2012338/h12v04_2012338_MCD_fsnow.tif': -1930181337,
        'modis/tiles/h12v04/2012338/h12v04_2012338_MCD_indices.tif': -521042334,
        'modis/tiles/h12v04/2012338/h12v04_2012338_MCD_snow.tif': 387672365,
        'modis/tiles/h12v04/2012338/h12v04_2012338_MOD-MYD_obstime.tif': -1046273330,
        'modis/tiles/h12v04/2012338/h12v04_2012338_MOD-MYD_temp.tif': 207791084,
        'modis/tiles/h12v04/2012338/h12v04_2012338_MOD_clouds.tif': 296967275,
    },
    'ignored': [ # contain site-specific data not applicable across users/configurations
        'gips-inv-db.sqlite3',
        'modis/tiles/h12v04/2012336/MCD43A2.A2012336.h12v04.006.2016112010833.hdf.index',
        'modis/tiles/h12v04/2012336/MCD43A4.A2012336.h12v04.006.2016112010833.hdf.index',
        'modis/tiles/h12v04/2012336/MOD10A1.A2012336.h12v04.005.2012339213007.hdf.index',
        'modis/tiles/h12v04/2012336/MOD11A1.A2012336.h12v04.006.2016131164025.hdf.index',
        'modis/tiles/h12v04/2012336/MYD10A1.A2012336.h12v04.005.2012340031954.hdf.index',
        'modis/tiles/h12v04/2012336/MYD11A1.A2012336.h12v04.006.2016132100529.hdf.index',
        'modis/tiles/h12v04/2012337/MCD43A2.A2012337.h12v04.006.2016112013509.hdf.index',
        'modis/tiles/h12v04/2012337/MCD43A4.A2012337.h12v04.006.2016112013509.hdf.index',
        'modis/tiles/h12v04/2012337/MOD09Q1.A2012337.h12v04.006.2015253024347.hdf.index',
        'modis/tiles/h12v04/2012337/MOD10A1.A2012337.h12v04.005.2012340033542.hdf.index',
        'modis/tiles/h12v04/2012337/MOD11A1.A2012337.h12v04.006.2016131185851.hdf.index',
        'modis/tiles/h12v04/2012337/MOD11A2.A2012337.h12v04.006.2016137164847.hdf.index',
        'modis/tiles/h12v04/2012337/MYD10A1.A2012337.h12v04.005.2012340112013.hdf.index',
        'modis/tiles/h12v04/2012337/MYD11A1.A2012337.h12v04.006.2016132101831.hdf.index',
        'modis/tiles/h12v04/2012338/MCD43A2.A2012338.h12v04.006.2016112020013.hdf.index',
        'modis/tiles/h12v04/2012338/MCD43A4.A2012338.h12v04.006.2016112020013.hdf.index',
        'modis/tiles/h12v04/2012338/MOD10A1.A2012338.h12v04.005.2012341091201.hdf.index',
        'modis/tiles/h12v04/2012338/MOD11A1.A2012338.h12v04.006.2016131192106.hdf.index',
        'modis/tiles/h12v04/2012338/MYD10A1.A2012338.h12v04.005.2012340142152.hdf.index',
        'modis/tiles/h12v04/2012338/MYD11A1.A2012338.h12v04.006.2016132112414.hdf.index',
    ],
}

# trailing whitespace and other junk characters are in current output
t_info = { 'stdout':  u"""\x1b[1mGIPS Data Repositories (v0.8.2)\x1b[0m
\x1b[1m
Modis Products v1.0.0\x1b[0m
\x1b[1m
Terra 8-day Products
\x1b[0m   ndvi8       Normalized Difference Vegetation Index: 250m
   temp8td     Surface temperature: 1km                
   temp8tn     Surface temperature: 1km                
\x1b[1m
Nadir BRDF-Adjusted 16-day Products
\x1b[0m   indices     Land indices                            
   quality     MCD Product Quality                     
\x1b[1m
Terra/Aqua Daily Products
\x1b[0m   fsnow       Fractional snow cover data              
   obstime     MODIS Terra/Aqua overpass time          
   snow        Snow and ice cover data                 
   temp        Surface temperature data                
\x1b[1m
Standard Products
\x1b[0m   clouds      Cloud Mask                              
   landcover   MCD Annual Land Cover                   
"""}

t_project = { 'created': {
    '0': None,
    '0/2012336_MCD_fsnow.tif': -1883071404,
    '0/2012336_MCD_indices.tif': 16775451,
    '0/2012336_MCD_quality.tif': -565438400,
    '0/2012336_MCD_snow.tif': -1824464052,
    '0/2012336_MOD-MYD_obstime.tif': 75984532,
    '0/2012336_MOD-MYD_temp.tif': -1641854393,
    '0/2012336_MOD_clouds.tif': -1957614367,
    '0/2012337_MCD_fsnow.tif': -856980949,
    '0/2012337_MCD_indices.tif': -943628281,
    '0/2012337_MCD_quality.tif': -567163117,
    '0/2012337_MCD_snow.tif': -1690607189,
    '0/2012337_MOD-MYD_obstime.tif': 1283853420,
    '0/2012337_MOD-MYD_temp.tif': 1191019971,
    '0/2012337_MOD_clouds.tif': -415873821,
    '0/2012337_MOD_ndvi8.tif': 486163471,
    '0/2012337_MOD_temp8td.tif': 1925913017,
    '0/2012337_MOD_temp8tn.tif': -1935381834,
    '0/2012338_MCD_fsnow.tif': -1017381876,
    '0/2012338_MCD_indices.tif': -2141895762,
    '0/2012338_MCD_quality.tif': -1145887707,
    '0/2012338_MCD_snow.tif': -319441628,
    '0/2012338_MOD-MYD_obstime.tif': 2120166986,
    '0/2012338_MOD-MYD_temp.tif': -1448563677,
    '0/2012338_MOD_clouds.tif': 1789735888,
}}

t_project_two_runs = t_project

t_project_no_warp = { 'created': {
    '0': None, # directory
    '0/2012336_MCD_fsnow.tif': -232655043,
    '0/2012336_MCD_indices.tif': 1748308005,
    '0/2012336_MCD_quality.tif': 1342857096,
    '0/2012336_MCD_snow.tif': 1704870455,
    '0/2012336_MOD-MYD_obstime.tif': -921877250,
    '0/2012336_MOD-MYD_temp.tif': -142389004,
    '0/2012336_MOD_clouds.tif': 792250507,
    '0/2012337_MCD_fsnow.tif': -118176399,
    '0/2012337_MCD_indices.tif': -558289948,
    '0/2012337_MCD_quality.tif': 1342857096,
    '0/2012337_MCD_snow.tif': -1562861219,
    '0/2012337_MOD-MYD_obstime.tif': -266130329,
    '0/2012337_MOD-MYD_temp.tif': 125915217,
    '0/2012337_MOD_clouds.tif': 1172608606,
    '0/2012337_MOD_ndvi8.tif': -996255186,
    '0/2012337_MOD_temp8td.tif': 1918564798,
    '0/2012337_MOD_temp8tn.tif': 1646469409,
    '0/2012338_MCD_fsnow.tif': -50404254,
    '0/2012338_MCD_indices.tif': 1361372393,
    '0/2012338_MCD_quality.tif': -429347844,
    '0/2012338_MCD_snow.tif': 415741551,
    '0/2012338_MOD-MYD_obstime.tif': -721900363,
    '0/2012338_MOD-MYD_temp.tif': -299932762,
    '0/2012338_MOD_clouds.tif': -1110899594,
}}

# TODO there should be something here but nothing is saved here during manual runs.
# See https://github.com/Applied-GeoSolutions/gips/issues/54
t_tiles = { 'created': {'h12v04': None}}

t_tiles_copy = { 'created': {
    'h12v04': None, # directory
    'h12v04/h12v04_2012336_MCD_fsnow.tif': 1284302156,
    'h12v04/h12v04_2012336_MCD_indices.tif': -2042919995,
    'h12v04/h12v04_2012336_MCD_quality.tif': 1121349116,
    'h12v04/h12v04_2012336_MCD_snow.tif': -2069225181,
    'h12v04/h12v04_2012336_MOD-MYD_obstime.tif': 808053323,
    'h12v04/h12v04_2012336_MOD-MYD_temp.tif': -1801734793,
    'h12v04/h12v04_2012336_MOD_clouds.tif': -221229092,
    'h12v04/h12v04_2012337_MCD_fsnow.tif': 1361214837,
    'h12v04/h12v04_2012337_MCD_indices.tif': -2147420472,
    'h12v04/h12v04_2012337_MCD_quality.tif': -1108654,
    'h12v04/h12v04_2012337_MCD_snow.tif': 1201721272,
    'h12v04/h12v04_2012337_MOD-MYD_obstime.tif': 673261584,
    'h12v04/h12v04_2012337_MOD-MYD_temp.tif': -1130082355,
    'h12v04/h12v04_2012337_MOD_clouds.tif': 1101505794,
    'h12v04/h12v04_2012337_MOD_ndvi8.tif': -181882164,
    'h12v04/h12v04_2012337_MOD_temp8td.tif': 868273975,
    'h12v04/h12v04_2012337_MOD_temp8tn.tif': 1173355207,
    'h12v04/h12v04_2012338_MCD_fsnow.tif': -647359984,
    'h12v04/h12v04_2012338_MCD_indices.tif': 1738864618,
    'h12v04/h12v04_2012338_MCD_quality.tif': 1379793143,
    'h12v04/h12v04_2012338_MCD_snow.tif': -1222056036,
    'h12v04/h12v04_2012338_MOD-MYD_obstime.tif': 411701973,
    'h12v04/h12v04_2012338_MOD-MYD_temp.tif': 1212632414,
    'h12v04/h12v04_2012338_MOD_clouds.tif': -2052728372,
}}

t_stats = { 'created': {
    'clouds_stats.txt': -142855826,
    'fsnow_stats.txt': 1649245444,
    'indices_stats.txt': -896610352,
    'ndvi8_stats.txt': -493982920,
    'obstime_stats.txt': -1865008021,
    'quality_stats.txt': -914294800,
    'snow_stats.txt': 239300424,
    'temp8td_stats.txt': 1298400117,
    'temp8tn_stats.txt': -312645073,
    'temp_stats.txt': 2580004,
}}
