/*!
 * \ingroup gip_landsat
 * \brief
 * Utility to preprocess landsat imagery
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
	bpo::options_description localopts("Volumetric water reflectance Options");
	localopts.add_options()
		("r", bpo::value<float>()->default_value(0.54), "Water-air reflection")
		("srband", bpo::value<unsigned char>()->default_value(5), "Spectral reflectrance band number")
		;
	gip::Options Opts(ac,av,localopts);
	if (Opts.Exit()) return 0;
	GDALAllRegister();
	// End boilerplate code **************************

	// Parameters
    float r = Opts.varmap["r"].as<float>();
    float p = Opts.varmap["p"].as<float>();
    float pp = Opts.varmap["pp"].as<float>();
    float n = Opts.varmap["n"].as<float>();
    float Q = Opts.varmap["Q"].as<float>();
    int srband = Opts.varmap["srband"].as<unsigned char>();
    float A = ((1-p)*(1-pp)) / (n*n);

	if (Opts.Verbose()) {
		cout << "Landsat Processing:" << endl;
        cout << "\tInput File: " << Opts.InputFile(0) << endl;
		cout << "\tOutput File: " << Opts.OutputFile(1) << endl;
	}

    // Test for right number of images
    GeoImage gimg(Opts.InputFile(0));
    GeoImage gimgout(Opts.OutputFile(0),gimg,gimg.NumBands()-2,GDT_Float32);
    // Bands to process
    vector<int> bands;
    for (int i=1;i<=gimg.NumBands()-1;i++) {
        if (i != srband) bands.push_back(i);
    }
    
    GeoRasterIO<float> img_sr( gimg[srband-1] );

    int outband(1);
    for (vector<int>::const_iterator iBand=bands.begin(); iBand!=bands.end(); iBand++) { 
        GeoRasterIO<float> img( gimg[*iBand-1] ); 
        GeoRasterIO<float> imgout( gimgout[outband-1] );
        CImg<float> cimg, cimgout, cimg_sr, cimgtmp;
        vector< box<point> > Chunks = img.Chunk();
	    if (Opts.Verbose())	cout << "\tProcessing band " << *iBand << " in " << Chunks.size() << " chunks:";
        int chunknum(0);
	    for (vector< box<point> >::const_iterator iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
	        if (Opts.Verbose()) cout << " " << (chunknum++)+1 << std::flush;
            cimg = img.Read(*iChunk) * 0.0001;
            cimgout = imgout.Read(*iChunk);
            cimg_sr = img_sr.Read(*iChunk) * 0.0001; 

            cimgtmp = (cimg - cimg_sr);
            cimgout = cimgtmp.div(A + r * Q * cimgtmp);

            imgout.Write(cimgout,*iChunk);
        }
	    if (Opts.Verbose()) cout << endl;
        outband++;
    }
	return 0;
}
