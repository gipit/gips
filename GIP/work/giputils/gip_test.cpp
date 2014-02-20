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
#include <gip/GeoRasterIO.h>

using namespace std;
using namespace gip;

int main (int ac, char* av[]) {
	// Begin boilerplate code ************************
	// Add program specific options
	//using namespace boost::program_options;
	//options_description localopts("Test Options");
	//localopts.add_options()
	//	("algtest", value<bool>()->implicit_value(true)->default_value(false), "Algorithm tests")
	//	("updatetest", value<bool>()->implicit_value(true)->default_value(false), "Update test");
	//gip::Options Opts(ac,av,localopts);
	//if (Opts.Exit()) return 0;
	GDALAllRegister();
	// End boilerplate code **************************

	//GeoImage test0("MCD43A4.A2006185.h09v04.005.2008138223410.tif", true);
	//test0[0].SetNoData(NULL);
	//cout << "test0: " << test0.Info() << endl;
    cout << "Opts" << endl;
    cout << Opts.InputFile() << endl;
	string newfname;
	cout << "GeoImage Open: " << Opts.InputFile() << endl;
	GeoImage test(Opts.InputFile(), true);
	cout << "test: " << test.Info() << endl;

    cout << "NoData test" << endl;
    cout << "NoData: " << test[0].NoData() << endl;
    test[0].ClearNoData();
    cout << "Cleared NoData...NoData = " << test[0].NoData() << endl;
    test[0].SetNoData(5);
    cout << "Set NoData...NoData = " << test[0].NoData() << endl;

    return 0;
	{
		cout << "GeoImage Copy constructor:" << endl;
		GeoImage copytest(test);
		cout << "copytest: " << copytest.Info(false) << endl;
	}
	cout << "copytest gone out of scope" << endl;
	cout << "test: " << test.Info(false);
	cout << "GeoImage Copy #2 (copytest)" << endl;
	GeoImage copytest(test);
	cout << "copytest: " << copytest.Info(false) << endl;

	// Create new file
	cout << "GeoImage Create:" << endl;
	GeoImage createtest("giptest_create", test, GDT_Float32, 1);
	newfname = createtest.Filename();
	createtest[0].SetDescription("createtest");
	//createtest[0].SetNoData(5);
	cout << "createtest: " << createtest.Info() << endl;

	// Assignment
	cout << "Assignment (test = createtest)" << endl;
	cout << "test before: " << test.Info(false);
	test = createtest;
	cout << "test after: " << test.Info(false);
	cout << "Assignment test #2 (test = copytest)" << endl;
	test = copytest;
	cout << "test: " << test.Info(false) << endl;

	cout << "Colors:" << endl;
	//test.SetColor(Color::Red, 1);
	try {
		cout << "Blue Band = " << test["Blue"].Info() << endl;
	} catch (...) {}
	try {
		cout << "Red Band:" << test["Red"].Info() << endl;
	} catch (...) {}

	//cout << "GeoImage Update:" << endl;
	// testing pam
	//test[0].SetDescription("writing band 1 of test");
	//test[0].UnsetNoData();

	bool update = Opts.varmap["updatetest"].as<bool>();
	if (update) {
		GeoImage updatetest(newfname, true);
		updatetest[0].SetDescription("updatetest");
		cout << "updatetest: " << updatetest.Info(true) << endl;
	}

	// Algorithms
    bool alg = Opts.varmap["algtest"].as<bool>();
    if (alg) {
	    cout << "CreateCopy (create new file copy)" << endl;
	    GeoImage newcopy = Copy(test,"giptest_createcopy.tif");
	    cout << "newcopy: " << newcopy.Info() << endl;

	    cout << "Indices test" << endl;
	    GeoImage indicestest = Indices(test,"giptest_indices.tif",true,true,true,true,true);
	    cout << "indicestest: " << indicestest.Info() << endl;

	    // Function tests
	    cout << "Functions (newraster = test[0] > 5.0)" << endl;
	    GeoRaster newraster = test[0] > 1000;
	    newcopy[0] = newraster;
	    cout << "newraster = " << newraster.Info() << endl;
	    cout << "newcopy = " << newcopy.Info() << endl;
	    GeoImage functiontest = Copy(newcopy,"giptest_function");
	    cout << "functiontest = " << functiontest.Info() << endl;
    }

	/*cout << "Chunking with Chunk Size = " << Options::ChunkSize() << endl;
	GeoProcess<double> Ptest(test[0]);
	vector< box<point> > Chunks = Ptest.Chunk();
	vector< box<point> >::const_iterator iChunk;
	cout << Chunks.size() << " chunks" << endl;
	int i(0);
	for (iChunk=Chunks.begin(); iChunk!=Chunks.end(); iChunk++) {
		cout << "Chunk " << i++ << endl;
		CImg<double> bandchunk( Ptest.Read(*iChunk) );
	}*/

	// Main algorithm
	//cout << "Input Files: " << endl; for (int i=0;i<Opts.InputFiles().size();i++) cout << "\t" << Opts.InputFiles(i) << endl;
	//cout << "Output Files: " << endl; for (int i=0;i<Opts.OutputFiles().size();i++) cout << "\t" << Opts.OutputFiles(i) << endl;
	//cout << "Bands: "; for (int i=0;i<Opts.Bands().size();i++) cout << " " << Opts.Bands(i); cout << endl;
	return 0;
}
