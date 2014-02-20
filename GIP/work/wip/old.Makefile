# Makefile for Geospatial Image Processing Library Utilities
#
# This library incorporates GDAL, for geospatial IO and reprojection, and CImg, an image processing library
#
# Author: Matt Hanson

TARGETS = Options GeoData GeoRaster GeoImage GeoAlgorithms

CXXFLAGS= -I/usr/include/gdal/ -I./ -O2
#LDFLAGS	= -L gip/ -L../gip/
LIBS = -lgdal -lboost_program_options -lboost_filesystem -lboost_system

#gip_blur: gip_blur.o
#	$(CC) $(LDFLAGS) $(LIBS) $@.o -o $@

all: $(TARGETS)
	make clean
	make $(TARGETS)
	

$(TARGETS) :
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(@).cpp -o $@ $(LIBS)

install:
	cp $(TARGETS) /usr/local/bin/

uninstall:
	for file in $(TARGETS);  do $(RM) /usr/local/bin/$$file; done

clean:
	$(RM) -rf $(TARGETS) $(LIB) *.o
