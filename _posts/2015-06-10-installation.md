---
layout: page
title: "Installation"
category: doc
date: 2015-06-10 14:41:53
order: 0
---


These installation notes are for Ubuntu 14.04.  While the easiest method is to install GIPS from PyPi, there are some system dependencies that need to be installed beforehand.

    # first install the UbuntuGIS repository if not already on system
    $ sudo apt-get install python-software-properties
    $ sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
    $ sudo apt-get update

    # install dependencies (for GIPPY)
    $ sudo apt-get install python-setuptools python-numpy python-gdal g++ libgdal1-dev gdal-bin libboost-dev swig2.0 swig



    # download and install GIPS from PyPi

    $ sudo pip install gips


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


#### Configuring Repositories

Installing GIPS only installs the Python packages and scripts, it does not configure any repositories. Before GIPS is used for the first time, gips_config is run to set up repository locations and preferences, as well as create any necessary directories.


Editing GIPS settings file

    1) Edit gips/settings.py to update with your environment
        - DATABASES: set connection settings and database if data tile vectors are stored in PostGIS
        - REPOS: set 'rootpath' for each dataset you want to enable
        - EMAIL: administrator email address

    2) Create repository directories
        - Make sure the directories in settings exist.
        - For each directory create 'tiles' and 'vectors' subdirectories
        - For each dataset put a shapefile 'tiles.shp' holding the tiles for the dataset







