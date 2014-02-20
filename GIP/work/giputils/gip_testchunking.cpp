/*!
 * \ingroup gip_test
 * \brief
 * Utility to geographically crop multiple images
 *
 * \author Matt Hanson <mhanson@appliedgeosolutions.com>
 * \date 2011-01-11
*/

#include <iostream>
#include <gip/gip.h>
#include <gip/GeoImage.h>
#include <gip/GeoAlgorithms.h>
#include <gip/GeoRasterIO.h>

using namespace std;
using namespace gip;

int main (int ac, char* av[]) {
	// Begin boilerplate code ************************
	// Add program specific options
	bpo::options_description localopts("Test Options");
	localopts.add_options()
		("ref", bpo::value<bool>()->implicit_value(true)->default_value(false), "Reference Counting test")
		("updatetest", bpo::value<bool>()->implicit_value(true)->default_value(false), "Update test");
	gip::Options Opts(ac,av,localopts);
	if (Opts.Exit()) return 0;
	GDALAllRegister();
	// End boilerplate code **************************

	//GeoImage test0("MCD43A4.A2006185.h09v04.005.2008138223410.tif", true);
	//test0[0].SetNoData(NULL);
	//cout << "test0: " << test0.Info() << endl;

	GeoImage test(Opts.InputFile(), true);
	cout << "test: " << test.Info() << endl;

	cout << "Chunking with Chunk Size = " << Options::ChunkSize() << endl;
	GeoRasterIO<double> Ptest(test[0]);
	vector< box<point> > Chunks = Ptest.Chunk();
	vector< box<point> >::const_iterator iChunk;
	cout << Chunks.size() << " chunks" << endl;
	int i(0);
	for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
		cout << "Chunk " << i++ << endl;
		CImg<double> bandchunk( Ptest.Read(*iChunk) );
	}

	return 0;
}
