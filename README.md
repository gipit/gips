
    **GIPIF**: Geospatial Image Processing and Inventory Framework

    Copyright (C) 2014 Matthew A Hanson

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program. If not, see <http://www.gnu.org/licenses/>

# GIPIF Installation

First step is to clone this repository to your local machine if you have not already done so.
Several packages are required for GIPIF. These notes are for Ubuntu systems.

    1) First install the UbuntuGIS Repository:
    $ sudo apt-get install python-software-properties
    $ sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
    $ sudo apt-get update

    2) Then install the required dependencies
    $ sudo apt-get install python-setuptools python-numpy python-gdal g++ libgdal1-dev gdal-bin libboost-dev swig2.0 swig
    These are required for Py6S
    $ sudo apt-get install python-scipy python-matplotlib

    3) Edit gipif/settings.py to update with your environment
        - DATABASES: set connection settings and database if data tile vectors are stored in PostGIS
        - REPOS: set 'rootpath' for each dataset you want to enable
        - EMAIL: administrator email address

    4) Create repository directories
        - Make sure the directories in settings exist.  
        - For each directory create 'tiles' and 'vectors' subdirectories
        - For each dataset put a shapefile 'tiles.shp' holding the tiles for the dataset

    5) Then install GIPIF
    $ ./setup.py install

Some datasets may require 6S to perform atmospheric correction. Follow the instructions here to install 6S:

    1) Download latest version of 6S (6SV1.1) to a working directory from http://6s.ltdri.org

    2) untar 6SV-1.1.tar
    $ tar xvf 6SV-1.1.tar

    3) Install fortran compiler
    $ sudo apt-get install gfortran

    3) Edit 6SV1.1/Makefile
        Replace this line:
            FC     = g77 $(FFLAGS)
        with this line:
            FC      = gfortran -std=legacy -ffixed-line-length-none $(FFLAGS)

To update GIPIF at a later date (will require root privileges)

    1) Get latest changes from git
    $ cd [gipif directory]
    $ git pull

    2) Run setup
    $ ./setup.py install

    If there is an error about a bad header, remove the previous GIPIF package
    $ ls /usr/local/lib/python2.7/dist-packages
        look for GIPIF package name
    $ rm -rf /usr/local/lib/python2.7/dist-packages/[GIPIF-package-name]
        Then go to step 2 again

## GIPIF Development Note

For developing GIPIF, it is recommended that you use a python virtual environment 
This allows multiple users on the same system to independently develop without 
collisions. If you are in a virtual environment (ve), install or develop will install
to the ve instead of the system
