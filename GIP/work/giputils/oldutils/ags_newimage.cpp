/*!
 * \ingroup ags_utilities
 * \brief
 * Create new GDAL image file
 *
 * \author Matt Hanson <mhanson@appliedgeosolutions.com>
 * \date 2011-01-11
*/

#include <agsOptions.h>
#include <agsGeoImage.h>
#include <agsGeoAlgorithm.h>

#include <iostream>
#include <string>
#include <vector>
//#include <map>

using namespace std;
using namespace ags;

int main (int ac, char* av[]) {
	// Setup program options and parse command line
	boost::program_options::options_description localopts("Create new file options");
	localopts.add_options()
		("numbands,n",bpo::value<int>()->default_value(1),"Number of bands in output")
		("xsize,x",bpo::value<int>()->default_value(1000),"XSize of output")
		("ysize,y",bpo::value<int>()->default_value(1000),"YSize of output")
		("datatype,d", bpo::value<string>()->default_value("Byte"),
               		"Output datatype: Byte, UInt16, Int16, UInt32, Int32, Float32, Float64");
	ags::agsOptions Opts(localopts);
	try {
		Opts.ParseCmdLine(ac,av);
	} catch (std::exception& err) {
		cout << "Error: " << err.what() << endl;
		return 1;
	}
	if (Opts.Exit()) return 0;
	// End parse command line
	GDALAllRegister();
	CPLSetConfigOption("GDAL_PAM_ENABLED",Opts.PAM() ? "TRUE" : "FALSE");
	typedef float T;

	// Options
	static map<string, GDALDataType> dtypes;
	dtypes["Byte"] = GDT_Byte; 
	dtypes["UInt16"] = GDT_UInt16; dtypes["Int16"] = GDT_Int16;
	dtypes["UInt32"] = GDT_UInt32; dtypes["Int32"] = GDT_Int32;
	dtypes["Float32"] = GDT_Float32; dtypes["Float64"] = GDT_Float64;
	string s = Opts.varmap["datatype"].as<string>();
	cout << "datatype = " << s << endl;
	GDALDataType t =  dtypes[Opts.varmap["datatype"].as<string>().c_str()];

	// Main algorithm
	if (Opts.InputFiles().empty()) {
		for (int i=0;i<Opts.OutputFiles().size();i++)
			GeoImage<T> Image(Opts.varmap["xsize"].as<int>(),Opts.varmap["ysize"].as<int>(),Opts.varmap["numbands"].as<int>(),
				Opts.OutputFiles(i), Opts.Format(), dtypes[Opts.varmap["datatype"].as<string>()]);
	} else {
		for (int i=0;i<Opts.OutputFiles().size();i++)
			GeoImage<T> Image( GeoImage<T>(Opts.InputFiles(i)),Opts.varmap["numbands"].as<int>(),
				dtypes[Opts.varmap["datatype"].as<string>()],
				Opts.OutputFiles(i),Opts.Format()
				);

	}
	
	return 0;
}
