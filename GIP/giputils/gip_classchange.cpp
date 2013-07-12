/*!
 * \ingroup gip_test
 * \brief
 * Utility to geographically crop multiple images
 *
 * \author Matt Hanson <mhanson@appliedgeosolutions.com>
 * \date 2011-01-11
*/

#include <iostream>
#include <gip/GeoImage.h>
#include <gip/GeoAlgorithms.h>
#include <gip/Options.h>

using namespace std;
using namespace gip;

int main (int ac, char* av[]) {
	// Begin boilerplate code ************************
	// Add program specific options
	using namespace boost::program_options;
	options_description localopts("Class Change Options");
	localopts.add_options()
		("change", value<string>(), "If threshold met change class M to N in image1")
		("gt", value<int>(), "Threshold compared against image2 (greater than)")
		("lt", value<int>(), "Threshold compared against image2 (less than)")
		("eq", value<int>(), "Threshold compared against image2 (equal to)")
		;
	gip::Options Opts(ac,av,localopts);
	if (Opts.Exit()) return 0;
	GDALAllRegister();
	// End boilerplate code **************************

	typedef boost::geometry::model::d2::point_xy<float> point;
	typedef boost::geometry::model::box<point> bbox;

	// Thresholds
    string op;
    if (Opts.varmap.count("gt")) {
        op = "gt";
    } else if (Opts.varmap.count("lt")) {
        op = "lt";
    } else if (Opts.varmap.count("eq")) {
        op = "eq";
    } else {
        cout << "No threshold given";
        //Opts.Help();
    }

    // Test for right number of images

    int threshold = Opts.varmap[op].as<int>();
    vector<unsigned int> classes = ParseToInts(Opts.varmap["change"].as<string>());

	if (Opts.Verbose()) {
		cout << "Class Change" << endl;
        cout << "\tChanging class " << classes[0] << " to " << classes[1] << " in " << Opts.InputFile(0) << endl;
		cout << "\tCriteria: " << Opts.InputFile(1) << " " << op << " " << threshold << endl;
		cout << "\tOutputting to " << Opts.OutputFile(0) << endl;
	}

    GeoImage gimg(Opts.InputFile(0));
    GeoRasterIO<unsigned char> img( gimg[0] );
    GeoRasterIO<float> opimg( GeoImage(Opts.InputFile(1))[0] );
    GeoRasterIO<unsigned char> imgout( GeoImage(Opts.OutputFile(0),gimg)[0] );

    CImg<unsigned char> cimg, cimgout;
    CImg<float> copimg;

    vector<bbox> Chunks = img.Chunk();
	if (Opts.Verbose())	cout << "\tProcessing " << Chunks.size() << " chunks:";
    int chunknum(0);
	for (vector< bbox >::const_iterator iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
	    if (Opts.Verbose()) cout << " " << (chunknum++)+1 << std::flush;
        cimg = img.Read(*iChunk);
        cimgout = cimg;
        copimg = opimg.Read(*iChunk);
        if (op == "gt") {
            copimg.threshold(threshold,false,true);
        } else if (op == "lt") {
            copimg.threshold(threshold)^=1;
        } else { // equal to
            copimg = (copimg.get_threshold(threshold,false,true) | (copimg.get_threshold(threshold)^=1))^=1;
        }

	    cimg_forXY(copimg,x,y) {
            if (copimg(x,y) == 1 && cimg(x,y) == classes[0]) {
                cimgout(x,y) = classes[1];
            }
        }

        imgout.Write(cimgout,*iChunk);
    }
	if (Opts.Verbose()) cout << endl;
	return 0;
}
