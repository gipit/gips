---
layout: page
title: "Installation"
category: doc
date: 2015-06-10 14:41:53
order: 0
---

First step is to clone this repository to your local machine if you have not already done so.
Several packages are required for GIPS. These notes are for Ubuntu systems.

    1) First install the UbuntuGIS Repository:
    $ sudo apt-get install python-software-properties
    $ sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
    $ sudo apt-get update

    2) Then install the required dependencies
    $ sudo apt-get install python-setuptools python-numpy python-gdal g++ libgdal1-dev gdal-bin libboost-dev swig2.0 swig
    These are required for Py6S
    $ sudo apt-get install python-scipy python-matplotlib

    3) Edit gips/settings.py to update with your environment
        - DATABASES: set connection settings and database if data tile vectors are stored in PostGIS
        - REPOS: set 'rootpath' for each dataset you want to enable
        - EMAIL: administrator email address

    4) Create repository directories
        - Make sure the directories in settings exist.
        - For each directory create 'tiles' and 'vectors' subdirectories
        - For each dataset put a shapefile 'tiles.shp' holding the tiles for the dataset

    5) Then install GIPS
    $ ./setup.py install

Some datasets may require 6S to perform atmospheric correction. Follow the instructions here to install 6S:

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


To update GIPS at a later date (will require root privileges)

    1) Get latest changes from git
    $ cd [gips directory]
    $ git pull

    2) Run setup
    $ ./setup.py install

    If there is an error about a bad header, remove the previous GIPS package
    $ ls /usr/local/lib/python2.7/dist-packages
        look for GIPS package name
    $ rm -rf /usr/local/lib/python2.7/dist-packages/[GIPS-package-name]
        Then go to step 2 again

