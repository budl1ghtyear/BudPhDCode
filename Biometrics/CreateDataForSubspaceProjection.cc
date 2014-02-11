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

SampleC<FilenameC> GetFilesFromList(FilenameC const &file, DirectoryC const &dir, StringC const &ext );
SampleC<SArray1dC<Point2dC> > ReadPointsFromVectorFile(FilenameC const &file);
int main(int nargs, char** nargv)
{	
	OptionC opt(nargs,nargv);
	DirectoryC featdir = opt.String("d","/vol/vssp/lip-tracking/VanillaFeatures/","Directory Containing Vanilla Features");
	FilenameC training = opt.String("f","/vol/vssp/lip-tracking/GMMExperiments/ProtocolFiles/C1traininglist","Filenames for use in projection");
	StringC ext = opt.String("e",".TSV","Feature file extension");
	FilenameC output = opt.String("o","/vol/vssp/lip-tracking/StatisticalShapeModelsOfLips/TrainingData/C1Raw","Output");
	opt.Check();
	
	SampleC<FilenameC> myfiles = GetFilesFromList(training,featdir,ext);
	cout<<"Files to be used are as follows: \n"<<myfiles<<endl;	
	SampleC<SArray1dC<Point2dC> > rawdata;
	for(SampleIterC<FilenameC> it(myfiles); it; it++)
		rawdata.Append(ReadPointsFromVectorFile(*it).Copy());
	//Now create output file
	ofstream wfile;
	wfile.open(output.chars());
	wfile << rawdata.Size() << "\n";
	for(SampleIterC<SArray1dC<Point2dC> > it(rawdata); it; it++)
	{	
		wfile << (*it) <<"\n";
	}
	wfile.close();
	return 0;
}

SampleC<FilenameC> GetFilesFromList(FilenameC const &file, DirectoryC const &dir, StringC const &ext )
{
	SampleC<FilenameC> files;
	ifstream list;
	list.open(file,ifstream::in);
	if(list.is_open())
	{
		while(! list.eof())
		{
			StringC filename, dummy;
			list >> dummy;
			if(dummy != "")
			{
				//now compile the filename
				filename = (StringC)dir + dummy + ext;
				files.Append((FilenameC)filename.Copy());
			}
		}
		list.close();
	}	
	return files;
}

SampleC<SArray1dC<Point2dC> > ReadPointsFromVectorFile(FilenameC const &file)
{
	SampleC<SArray1dC<Point2dC> > points;
	ifstream myfile;
	myfile.open(file,ifstream::in);
	if(myfile.is_open())
	{
		UIntT numel(0),vectsize(0);
		myfile >> numel;
		myfile >> vectsize;
		if(vectsize != 22)
			cout<<"This file has vectors less than 22 in length: "<<file<<endl;
		for(UIntT i = 0; i < numel; i++)
		{
			VectorC vect(vectsize); vect.Fill(0.0);
			for(UIntT j = 0; j < vectsize; j++)
			{
				myfile >> vect[j];
			}
			points.Append(VectorToPoints(vect).Copy());
		}
		myfile.close();
	}
	return points;
}

