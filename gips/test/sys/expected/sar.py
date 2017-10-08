"""Known-good outcomes for tests, mostly stdout and created files."""

t_info = {
    'stdout':  u"""\x1b[1mGIPS Data Repositories (v0.8.2)\x1b[0m
\x1b[1m
SAR Products v0.9.0\x1b[0m
\x1b[1m
Standard Products
\x1b[0m   date        Day of year array                       
   linci       Incident angles                         
   sign        Sigma nought (radar backscatter coefficient)
"""
}


t_inventory = {
    'stdout': u"""\x1b[1mGIPS Data Inventory (v0.8.2)\x1b[0m
Retrieving inventory for site tiles

\x1b[1mAsset Holdings for tile N07E099\x1b[0m
\x1b[1m\x1b[4m    DATE      alos1     alos2    Product  \x1b[0m
\x1b[1m2015        
\x1b[0m    101                 100.0%   


1 files on 1 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mAFBD: PALSAR FineBeam Dual Polarization\x1b[0m
\x1b[31mAFBS: PALSAR FineBeam Single Polarization\x1b[0m
\x1b[32mAWB1: PALSAR WideBeam (ScanSAR Short Mode)\x1b[0m
\x1b[34mAWBD: PALSAR-2 WideBeam (ScanSAR w/Ortho Slope Correction)\x1b[0m
\x1b[33mJFBS: JERS-1 FineBeam Single Polarization\x1b[0m
Retrieving inventory for site tiles

\x1b[1mAsset Holdings for tile N19E100\x1b[0m
\x1b[1m\x1b[4m    DATE      alos1     alos2    Product  \x1b[0m
\x1b[1m2010        
\x1b[0m    182       100.0%             


1 files on 1 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mAFBD: PALSAR FineBeam Dual Polarization\x1b[0m
\x1b[31mAFBS: PALSAR FineBeam Single Polarization\x1b[0m
\x1b[32mAWB1: PALSAR WideBeam (ScanSAR Short Mode)\x1b[0m
\x1b[34mAWBD: PALSAR-2 WideBeam (ScanSAR w/Ortho Slope Correction)\x1b[0m
\x1b[33mJFBS: JERS-1 FineBeam Single Polarization\x1b[0m
Retrieving inventory for site tiles

\x1b[1mAsset Holdings for tile N00E099\x1b[0m
\x1b[1m\x1b[4m    DATE      alos1     alos2    Product  \x1b[0m
\x1b[1m2009        
\x1b[0m    041       100.0%             


1 files on 1 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mAFBD: PALSAR FineBeam Dual Polarization\x1b[0m
\x1b[31mAFBS: PALSAR FineBeam Single Polarization\x1b[0m
\x1b[32mAWB1: PALSAR WideBeam (ScanSAR Short Mode)\x1b[0m
\x1b[34mAWBD: PALSAR-2 WideBeam (ScanSAR w/Ortho Slope Correction)\x1b[0m
\x1b[33mJFBS: JERS-1 FineBeam Single Polarization\x1b[0m
""",
}

otro = u"""\x1b[1mGIPS Data Inventory (v0.8.2)\x1b[0m
Retrieving inventory for site tiles

\x1b[1mAsset Holdings for tile N07E099\x1b[0m
\x1b[1m\x1b[4m    DATE      alos1     alos2    Product  \x1b[0m
\x1b[1m2015        
\x1b[0m    101                 100.0%   


1 files on 1 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mAFBD: PALSAR FineBeam Dual Polarization\x1b[0m
\x1b[31mAFBS: PALSAR FineBeam Single Polarization\x1b[0m
\x1b[32mAWB1: PALSAR WideBeam (ScanSAR Short Mode)\x1b[0m
\x1b[34mAWBD: PALSAR-2 WideBeam (ScanSAR w/Ortho Slope Correction)\x1b[0m
\x1b[33mJFBS: JERS-1 FineBeam Single Polarization\x1b[0m
"""


t_process = {
    'updated': {
        'sar/stage': None,
        'sar/tiles/N00E099/2009041': None,
        'sar/tiles/N07E099/2015101': None,
        'sar/tiles/N19E100/2010182': None,
    },
    'created': {
        'sar/tiles/N00E099/2009041/KC_017-C25N00E099WB1ORSA1.tar.gz.index': -741988802,
        'sar/tiles/N00E099/2009041/N00E099_2009041_AWB1_date.tif': 307502521,
        'sar/tiles/N00E099/2009041/N00E099_2009041_AWB1_linci.tif': -399377879,
        'sar/tiles/N00E099/2009041/N00E099_2009041_AWB1_sign.tif': 1287316543,
        'sar/tiles/N00E099/2009041/N00E099_K25_date': 717902658,
        'sar/tiles/N00E099/2009041/N00E099_K25_date.hdr': 886355476,
        'sar/tiles/N07E099/2015101/KC_999-C019DRN07E099WBDORSA1.tar.gz.index': 541661026,
        'sar/tiles/N07E099/2015101/N07E099_2015101_AWBD_date.tif': 1151011840,
        'sar/tiles/N07E099/2015101/N07E099_2015101_AWBD_linci.tif': -2017292803,
        'sar/tiles/N07E099/2015101/N07E099_2015101_AWBD_sign.tif': 515521946,
        'sar/tiles/N07E099/2015101/N07E099_W02DC019DR_date': 1689631779,
        'sar/tiles/N07E099/2015101/N07E099_W02DC019DR_date.hdr': 1634871345,
        'sar/tiles/N19E100/2010182/KC_017-Y10N19E100FBDORSA1.tar.gz.index': 15900329,
        'sar/tiles/N19E100/2010182/N19E100_10_date': 330702546,
        'sar/tiles/N19E100/2010182/N19E100_10_date.hdr': -1922938697,
        'sar/tiles/N19E100/2010182/N19E100_2010182_AFBD_date.tif': 238211992,
        'sar/tiles/N19E100/2010182/N19E100_2010182_AFBD_linci.tif': 1373550692,
        'sar/tiles/N19E100/2010182/N19E100_2010182_AFBD_sign.tif': 832541493
    },
    '_symlinks': {
    },
    'ignored': [
        'gips-inv-db.sqlite3',
    ],
    '_inv_stdout': """\x1b[1mGIPS Data Inventory (v0.8.2)\x1b[0m
Retrieving inventory for site tiles

\x1b[1mAsset Holdings for tile N07E099\x1b[0m
\x1b[1m\x1b[4m    DATE      alos1     alos2    Product  \x1b[0m
\x1b[1m2015        
\x1b[0m    101                 100.0%     \x1b[34mdate\x1b[0m  \x1b[34mlinci\x1b[0m  \x1b[34msign\x1b[0m


1 files on 1 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mAFBD: PALSAR FineBeam Dual Polarization\x1b[0m
\x1b[31mAFBS: PALSAR FineBeam Single Polarization\x1b[0m
\x1b[32mAWB1: PALSAR WideBeam (ScanSAR Short Mode)\x1b[0m
\x1b[34mAWBD: PALSAR-2 WideBeam (ScanSAR w/Ortho Slope Correction)\x1b[0m
\x1b[33mJFBS: JERS-1 FineBeam Single Polarization\x1b[0m
Retrieving inventory for site tiles

\x1b[1mAsset Holdings for tile N19E100\x1b[0m
\x1b[1m\x1b[4m    DATE      alos1     alos2    Product  \x1b[0m
\x1b[1m2010        
\x1b[0m    182       100.0%               \x1b[35mdate\x1b[0m  \x1b[35mlinci\x1b[0m  \x1b[35msign\x1b[0m


1 files on 1 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mAFBD: PALSAR FineBeam Dual Polarization\x1b[0m
\x1b[31mAFBS: PALSAR FineBeam Single Polarization\x1b[0m
\x1b[32mAWB1: PALSAR WideBeam (ScanSAR Short Mode)\x1b[0m
\x1b[34mAWBD: PALSAR-2 WideBeam (ScanSAR w/Ortho Slope Correction)\x1b[0m
\x1b[33mJFBS: JERS-1 FineBeam Single Polarization\x1b[0m
Retrieving inventory for site tiles

\x1b[1mAsset Holdings for tile N00E099\x1b[0m
\x1b[1m\x1b[4m    DATE      alos1     alos2    Product  \x1b[0m
\x1b[1m2009        
\x1b[0m    041       100.0%               \x1b[32mdate\x1b[0m  \x1b[32mlinci\x1b[0m  \x1b[32msign\x1b[0m


1 files on 1 dates
\x1b[1m
SENSORS\x1b[0m
\x1b[35mAFBD: PALSAR FineBeam Dual Polarization\x1b[0m
\x1b[31mAFBS: PALSAR FineBeam Single Polarization\x1b[0m
\x1b[32mAWB1: PALSAR WideBeam (ScanSAR Short Mode)\x1b[0m
\x1b[34mAWBD: PALSAR-2 WideBeam (ScanSAR w/Ortho Slope Correction)\x1b[0m
\x1b[33mJFBS: JERS-1 FineBeam Single Polarization\x1b[0m
""",
}

t_project = {
    'created': {
    }
}

t_project_no_warp = {
    'created': {'0': None,
    }
}


t_tiles_copy = {
    'created': {
    }
}

t_tiles_copy = {
    'created': {
    }
}

t_stats = { 
    'created': {
    }
}

