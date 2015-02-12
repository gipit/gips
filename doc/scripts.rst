*******************************************************************************
GIPS Scripts
*******************************************************************************

:Author: Matthew Hanson
:Contact: matt.a.hanson@gmail.com

While GIPS can be used as a library in existing applications, it is mainly used as a series of command line scripts that allow users to query, manage, download and process raster imagery.

Most scripts require a dataset name, specified like a command (e.g., gips_info Landsat). These are indicated as below with *dataname* following the script name. The dataset name is case sensitive, and a complete list of available datasets will be printed out when calling any script without a dataset name.


``gips_info`` *dataset*
------------------------------------------------------------------------------
*gips_info* will output a list of products available for the given dataset. 

``gips_inventory`` *dataset*
------------------------------------------------------------------------------
*gips_inventory* provides the basic functionality for creating a data inventory and printing it. An inventory is a query of what is currently available in a data repository (i.e., a dataset). An inventory is also implicitly created by most other scripts, and thus the options available to *gips_inventory* are all available for other commands as well.
