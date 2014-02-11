//This program aims to unbundle the spl biometrics lip feature data created by URS into individual TSV files

//INCLUDE FILES
#include <iostream>
#include <fstream>
#include 	"Ravl/Option.hh"
#include 	"Ravl/OS/Filename.hh"
#include 	"Ravl/OS/Directory.hh"
#include 	"Ravl/DLIter.hh"
#include 	"Ravl/DList.hh"
#include 	"Ravl/SArray1d.hh"
#include	"Ravl/SArray1dIter.hh"
#include 	"Ravl/PatternRec/Sample.hh"
#include 	"Ravl/PatternRec/SampleIter.hh"
#include 	"Ravl/Vector.hh"
#include 	"BANCAGT.hh"
#include 	"UtilityFunctions.hh"
#include 	"MBspline.hh"
using namespace RavlN;
using namespace std;
//MAIN PROGRAM

int main(int nargs, char** nargv)
{	
	OptionC opt(nargs,nargv);
	FilenameC training = opt.String("f","/vol/vssp/lip-tracking/GMMExperiments/ProtocolFiles/C1traininglist","Filenames for use in projection");
	FilenameC output = opt.String("o","/vol/vssp/lip-tracking/StatisticalShapeModelsOfLips/TrainingData/TSVC1Raw","Output");
	opt.Check();
	
	SampleC<SArray1dC<Point2dC> > data;
	IStreamC file(training);
	file >> data;
	cout<<"Size of sample is "<<data.Size()<<endl;
	ofstream wfile;
	wfile.open(output.chars());
	wfile << data.Size() << "\t" << data[0].Size()*data[0][0].Size() << "\n";

	for(SampleIterC<SArray1dC<Point2dC> > it(data); it; it++)
	{
		VectorC vect = PointsToVector(*it).Copy();
		for(UIntT i = 0; i < vect.Size(); i++)
			wfile<<vect[i]<<"\t";
		wfile<<"\n";
	}
	wfile.close();
	return 0;
}

