/*!
 * \ingroup ags_utilities
 * \brief
 * Utility to geographically crop multiple images
 *
 * \author Matt Hanson <mhanson@appliedgeosolutions.com>
 * \date 2011-01-11
*/

// STL
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include <agsOptions.h>
#include <agsGeoImage.h>

using namespace std;
using namespace ags;

int main (int ac, char* av[]) {
	// Setup program options and parse command line
	boost::program_options::options_description localopts("Test Options");
	localopts.add_options()
		("test,t", bpo::value<int>(), "Print values only");
	ags::agsOptions Opts(localopts);
	try {
		Opts.ParseCmdLine(ac,av);
	} catch (std::exception& err) {
		cout << "Error: " << err.what() << endl;
		return 1;
	}
	// End parse command line
	cout << "Input Files: "; for (int i=0;i<Opts.InputFiles().size();i++) cout << Opts.InputFiles(i) << " "; cout << endl;

	GDALAllRegister();

	// Main algorithm
	GeoImage<unsigned char> ImageOut;
	unsigned char w1, w2, w3, w4;
	for(int iFiles=0;iFiles<Opts.InputFiles().size();iFiles++) {
		// Read in this input image
		cout << "read input" << endl;
		GeoImage<unsigned char> ImageIn(Opts.InputFiles(iFiles));
		
		if (iFiles == 0) {
			ImageOut = GeoImage<unsigned char>(ImageIn,Opts.InputFiles().size(),Opts.OutputFile(),Opts.Format());
		}

		vector< bgeo::box<point_2d> > Chunks = ImageIn.Chunk();
		CImg<unsigned char> imgout, w1img, w2img, w3img, w4img;
		cout << "chunk loop" << endl;
		for (int chunk=0;chunk<Chunks.size(); chunk++) {
                        w1img = ImageIn[0].Read(Chunks[chunk]);
                        w2img = ImageIn[1].Read(Chunks[chunk]);
                        w3img = ImageIn[2].Read(Chunks[chunk]);
                        w4img = ImageIn[3].Read(Chunks[chunk]);
                        imgout = ImageOut[iFiles].Read(Chunks[chunk]);

			cimg_forXY(imgout,x,y) {
				w1=w1img(x,y);
				w2=w2img(x,y);
				w3=w3img(x,y);
				w4=w4img(x,y);
				//if ( ((w1==w2) && (w1==1)) || ((w2==w3) && (w2==1)) || ((w3==w4) && (w3==1)) ) imgout(x,y) = 1;
				//if ((w1+w2+w3+w4 > 2) && (w2 != 0) && (w3 != 0)) imgout(x,y) = 1;
				if ((w1+w2+w3+w4 > 3)) imgout(x,y) = 1;
			}

                        // Writing Result
                        ImageOut[iFiles].Write(imgout,Chunks[chunk]);
                }
		//ImageOut[iFiles].SetDescription( ImageIn.Basename().c_str() );
		//ImageOut[iFiles].GetGDALDataset()->FlushCache();
	}
}
