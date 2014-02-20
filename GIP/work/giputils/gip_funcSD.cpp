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
#include <gip/GeoImageIO.h>
#include <gip/GeoAlgorithms.h>

using namespace std;
using namespace gip;

int main (int ac, char* av[]) {
	// Begin boilerplate code ************************
	// Add program specific options
	bpo::options_description localopts("Apply function");
	localopts.add_options()
		("func", bpo::value<string>(), "Function to apply to bands")
		;
	gip::Options Opts(ac,av,localopts);
	if (Opts.Exit()) return 0;
	GDALAllRegister();
	// End boilerplate code **************************

	// Parameters
    //string func = Opts.varmap["func"].as<string>();

	if (Opts.Verbose()) {
		cout << "Band math function" << endl;
        cout << "\tInput File: " << Opts.InputFile(0) << endl;
		cout << "\tOutput File: " << Opts.OutputFile(0) << endl;
	}

    // Parse formula
    vector<int> bands;

    // Test for right number of images
    GeoImage gimg(Opts.InputFile(0));
    GeoImage gimgout(Opts.OutputFile(0),gimg,1,GDT_Float32);

    GeoRasterIO<float> img1( gimg[0] );
    GeoRasterIO<float> img3( gimg[2] );
    GeoRasterIO<float> imgout(gimgout);
    CImg<float> cimg1, cimg3, cimgout;

    vector< box<point> > Chunks = img1.Chunk();
    int chunknum(0);
    for (vector< box<point> >::const_iterator iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
        if (Opts.Verbose()) cout << " " << (chunknum++)+1 << std::flush;
        cimg1 = img1.Read(*iChunk); 
        cimg3 = img3.Read(*iChunk);
        cimgout = imgout.Read(*iChunk); 
       
        // Algorithm 
        cimgout = (6.155*cimg1) + (-14.357*cimg3) + 1.8075;
        cimgout = cimgout.mul(cimgout);

        imgout.Write(cimgout,*iChunk); 
    }
	return 0;
}
