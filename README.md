
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

Several packages are required for GIPIF

1) First install the UbuntuGIS Repository:
$ sudo apt-get install python-software-properties
$ sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
$ sudo apt-get update

2) Then install the required dependencies
$ sudo apt-get install python-setuptools python-numpy python-gdal g++ libgdal1-dev gdal-bin libboost-dev-all swig2.0 swig

3) Then install GIPIF
$ ./setup.py install


## GIPIF Development Note

For developing GIPIF, it is recommended that you use a python virtual environment 
This allows multiple users on the same system to independently develop without 
collisions. If you are in a virtual environment (ve), install or develop will install
to the ve instead of the system
