/*!
 * \ingroup gip_test
 * \brief
 * Utility to geographically crop multiple images
 *
 * \author Matt Hanson <mhanson@appliedgeosolutions.com>
 * \date 2011-01-11
*/

#include <iostream>
#include <gip/GeoImageIO.h>
#include <gip/GeoAlgorithms.h>
#include <gip/Options.h>

using namespace std;
using namespace gip;

int main (int ac, char* av[]) {
	// Begin boilerplate code ************************
	// Add program specific options
	using namespace boost::program_options;
	using boost::geometry::model::box;
	options_description localopts("Rice Detect Options");
	localopts.add_options()
		("maxcrops", value<int>()->default_value(0), "Maximum number of crops to record (i.e. peaks)")
		("thlow", value<double>()->default_value(0), "Low Threshold (<=)")
		("thhigh", value<double>()->default_value(0), "High Threshold (>=)")
		("dthlow", value<int>()->default_value(90), "Low Day Threshold (<=)")
		("dthhigh", value<int>()->default_value(120), "High Day Threshold (>=)")
		;
	gip::Options Opts(ac,av,localopts);
	if (Opts.Exit()) return 0;
	GDALAllRegister();
	// End boilerplate code **************************

	typedef float T;

	int maxcrops = Opts.varmap["maxcrops"].as<int>();
	int numbands = 2 * maxcrops + 3;
	//int numbands(maxcrops +2);
	// Thresholds
	double th0( Opts.varmap["thlow"].as<double>() );
	double th1( Opts.varmap["thhigh"].as<double>() );
	int dth0( Opts.varmap["dthlow"].as<int>() );
    int dth1( Opts.varmap["dthhigh"].as<int>() );

	if (Opts.Verbose()) {
		cout << "AGS Rice Detection" << endl;
		cout << "\tLow Threshold = " << th0 << endl;
		cout << "\tHigh Threshold = " << th1 << endl;
		cout << "\tDate Range = " << dth0 << " - " << dth1 << endl;
		cout << "\tMax Crops = " << maxcrops << endl;
		cout << "\t" << Opts.InputFiles().size() << " files to process" << endl;
	}
	for (unsigned int f=0;f<Opts.InputFiles().size();f++) {
		// for loop here to go through all files
		GeoImageIO<T> Image(GeoImage(Opts.InputFile(f)));
		GeoImageIO<unsigned char> ImageOut(GeoImage(Opts.OutputFile(f), Image, GDT_Byte, numbands));

		// Get timestamp, in days after first image
		vector<string> bandnames = Image.BandNames();
		vector<int> doy;
		for (unsigned int i=0; i<bandnames.size(); i++) {
			//doy.push_back( atoi(bandnames[i].c_str()) - atoi(bandnames[0].c_str()) );
			doy.push_back( atoi(bandnames[i].c_str()) ); // - atoi(bandnames[0].c_str()) );
		}

		if (Opts.Verbose()) {
			cout << Opts.InputFile(f) << " -> " << Opts.OutputFile(f) << endl;
			cout << "\tTemporal Bands (DOY): ";
			for (unsigned int i=0;i<doy.size();i++) cout << " " << doy[i];
			cout << endl;
		}

		point p1, p2;
		int width, height;

		//CImg<unsigned int> Clow( iChunk->width(), iChunk->height(), 1, 1, 0 );
		vector< box<point> > Chunks = Image.Chunk();
		int chunknum(0);
		int peaknum(0);
		if (Opts.Verbose())	cout << "\tProcessing " << Chunks.size() << " chunks:";
		for (vector< box<point> >::const_iterator iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
			if (Opts.Verbose()) cout << " " << (chunknum++)+1 << std::flush;
			CImgList<T> CImage( Image.ReadAsList(*iChunk) );
			CImgList<T> CImageOut( ImageOut.ReadAsList(*iChunk) );

			p1 = iChunk->min_corner();
			p2 = iChunk->max_corner();
			width = p2.x()-p1.x()+1;
			height = p2.y()-p1.y()+1;

			// Reset running DOY to all zero
			CImg<int> DOY( width, height, 1, 1, 0 );

			// Seed with first band
			CImg<int> Clow, Chigh, matched;
			CImg<int> ClowEver( CImage[0].get_threshold(th0)^1 );
			CImg<int> ChighEver( CImage[0].get_threshold(th1,false,true) );
			CImageOut[1] = ClowEver;
			CImageOut[2] = ChighEver;
			for (unsigned int b=1;b<Image.NumBands();b++) {
				// Get low and high thresholds for this band
				Clow = CImage[b].get_threshold(th0);
				Chigh = CImage[b].get_threshold(th1,false,true);

				// Temporal processing
				// Where <= low threshold, set DOY to 0
				DOY.mul(Clow);
				// Where > low threshold, add to previous DOY
				DOY += Clow * doy[b];

				// Update flags if it was Ever low/high
				ClowEver |= (Clow^=1); // Clow now represents <= th0
				ChighEver |= Chigh;

				// Add to number of troughs/peaks
				CImageOut[1] += Clow;
				CImageOut[2] += Chigh;

				// More temporal processing
				// If low threshold was never met, change DOY to zero
				DOY.mul(ClowEver^1);
				Chigh = Chigh & DOY.get_threshold(dth0,true) & (DOY.get_threshold(dth1)^1);
				matched += Chigh;
				// Reset if high date has passed
				DOY.mul(DOY.get_threshold(dth1));
				// Loop through
				if (maxcrops != 0) {
					cimg_forXY(CImageOut[0],x,y) {
						if (Clow(x,y)) {
							peaknum = CImageOut[1](x,y);
							if (peaknum > 0 && peaknum <= maxcrops) {
								CImageOut[(peaknum-1)*2+3](x,y) = doy[b];
							}
						}
						if (Chigh(x,y)) {
							peaknum = CImageOut[2](x,y);
							if (peaknum > 0 && peaknum <= maxcrops) {
								CImageOut[(peaknum-1)*2+4](x,y) = doy[b];
								// Day of Year
								CImageOut[(peaknum-1)*2+3](x,y) = DOY(x,y); // * Chigh(x,y);
								// Yield (value)
								CImageOut[(peaknum-1)*2+4](x,y) = CImage[b](x,y);
							}
						}

					}
				}
			}
			CImageOut[0] = ClowEver & ChighEver;
			CImg<T> tmp = CImageOut[0].get_threshold(0);
			for (unsigned int c=0;c<CImageOut.size();c++) {
				CImageOut[c].mul(tmp);
			}

			// Write out images
			ImageOut.Write(CImageOut.get_append('c'), *iChunk, true);
		}
		if (Opts.Verbose()) cout << endl;
	}
	return 0;
}
