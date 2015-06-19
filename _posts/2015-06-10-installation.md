---
layout: page
title: "Installation"
category: doc
date: 2015-06-10 14:41:53
order: 0
---


These installation notes are for Ubuntu 14.04.  While the easiest method is to install GIPS from PyPi using pip, there are some system dependencies that need to be installed beforehand.

    # install UbuntuGIS repository
    $ sudo apt-get install python-software-properties
    $ sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
    $ sudo apt-get update

    # install dependencies
    $ sudo apt-get install python-setuptools python-numpy python-gdal g++ libgdal1-dev gdal-bin libboost-dev swig2.0 swig

    # install GIPS from PyPi
    $ sudo pip install gips

This installs the GIPS package and scripts on the system, but data repositories must first be set up in order to use GIPS.


### Configuring Repositories
Installing GIPS only installs the Python packages and scripts, it does not configure any repositories. Before GIPS is used for the first time, gips_config is run to set up repository locations and preferences, as well as create any necessary directories.


#### Environment settings
~~~
# Used for anonymous FTP
EMAIL = '$EMAIL'


# Site files and data tiles vectors can be retrieved from a database
DATABASES = {
#    'tiles': {
#        'NAME': '',
#        'USER': '',
#        'PASSWORD': '',
#        'HOST': '',
#        'PORT': '5432',
#    }
}


REPOS = {
    'aod': {
        'repository': '$TLD/aod',
    },
    'landsat': {
        'repository': '$TLD/landsat',
        # Landsat specific settings
        '6S': False,            # atm correction for VIS/NIR/SWIR bands
        'MODTRAN': False,       # atm correction for LWIR
        'extract': False,       # extract files from tar.gz before processing instead of direct access
    },
    'modis': {
        'repository': '$TLD/modis',
    },
    # these drivers tend to more specialized and experimental so turned off by default
    #'cdl': {
    #    'repository': '$TLD/cdl',
    #},
    #'sar': {
    #    'repository': '$TLD/sar',
    #},
    #'sarannual': {
    #    'repository': '$TLD/sarannual',
    #},
    #'merra': {
    #    'repository': '$TLD/Merra',
    #},
    #'daymet': {
    #    'repository': '$TLD/daymet',
    #},
}


"""
# to add repository add new key to the REPOS dictionary
    'dataname': {
        # path to driver directory location (default to gips/data/dataname/ if not given)
        'driver': '',
        # path to top level directory of data
        'repository': '',
        # override location of tiles vector (default to gips/data/dataname/tiles.shp)
       'tiles': '',
        #'tiles': 'mydatabase:mydatatype_tiles',        # database format
        #'tiles': '~/randomdir/dataname_tiles.shp'      # file format
    }
"""
~~~


#### User settings
~~~
# this block will try to read in an environment (e.g., system-level or virtual)
# settings file so that individual settings can be changed
try:
    import gips.settings
    execfile(gips.settings.__file__.rstrip('c'))
except:
    raise Exception('There are no environment-level GIPS settings!')


# change email used for FTP
#EMAIL = ''


# to add in a database add a new key to the DATABASES dictionary
#DATABASES['mydatabase'] = {
#        'NAME': '',
#        'USER': '',
#        'PASSWORD': '',
#        'HOST': '',
#        'PORT': '5432',
#}


"""
# to add repository add new key to the REPOS dictionary
REPOS['dataname'] = {
    # path to driver location (default to gips/data/dataname)
    'driver': '',
    # path to top level directory of data
    'repopath': '',
    # override location of tiles vector (defaults to gips/data/dataname/tiles.shp)
   'tiles': '',
    #'tiles': 'mydatabase:mydatatype_tiles',        # database format
    #'tiles': '~/randomdir/dataname_tiles.shp'      # file format

   # 'attribute name holding tileid in tiles vector'
   'tileid_attribute': '',
}
"""
~~~

    1) Edit gips/settings.py to update with your environment
        - DATABASES: set connection settings and database if data tile vectors are stored in PostGIS
        - REPOS: uncomment dataset to enable
        - EMAIL: administrator email address (currently used for FTP)

If the file is edited, the gips_config command should be run again, as this will make sure all the necessary directories exist for each enabled repository.


#### Atmospheric Correction Dependencies
Currently the Landsat driver requires 6S to perform atmospheric correction in the visible bands. Follow the instructions here to install 6S:

    1) Download latest version of 6S (6SV1.1) to a working directory from http://6s.ltdri.org

    2) untar 6SV-1.1.tar
    $ tar xvf 6SV-1.1.tar

    3) Install fortran compiler and make
    $ sudo apt-get install gfortran make

    3) Edit 6SV1.1/Makefile
        Replace this line:
            FC     = g77 $(FFLAGS)
        with this line:
            FC      = gfortran -std=legacy -ffixed-line-length-none $(FFLAGS)

    4) Compile 6S
    $ cd 6SV1.1
    $ make

    These are required for Py6S
    $ sudo apt-get install python-scipy python-matplotlib