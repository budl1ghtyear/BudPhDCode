//This program aims to produce a set of PCA projected feature vectors given some TSV vanilla features of some subjects

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
#include 	"ICAProjection.hh"
using namespace RavlN;
using namespace std;
//MAIN PROGRAM

void WriteVectorsToTSVFile(SampleC<VectorC> const &pts, FilenameC const &file, FilenameC const &origfile);
SampleC<SArray1dC<Point2dC> > ReadPointsFromVectorFile(FilenameC const &file);
int main(int nargs, char** nargv)
{	
	OptionC opt(nargs,nargv);
	DirectoryC icadir = opt.String("p","/vol/vssp/lip-tracking/StatisticalShapeModelsOfLips/C1TrainingModels/ICA/","Directory containing pca data");
	DirectoryC vanilla = opt.String("d","/vol/vssp/lip-tracking/VanillaFeatures/","Directory containing vanilla features");
	StringC ext = opt.String("e",".TSV","Extension of the vanilla features");
	DirectoryC output = opt.String("o","/vol/vssp/lip-tracking/PCAShapeFeatures/","Output directory in which to store the file");
	UIntT numcomp = opt.Int("n",10,"Number of eigencomponents to preserve, default = 10");
	opt.Check();
	
	//LOAD PCA DATA
	ICAProjectionC pca((StringC)icadir);
	StringC extension = "*" + ext;
	DListC<StringC> files = vanilla.FiltList(extension);
	for(DLIterC<StringC> it(files); it; it++)
	{
		StringC infile = (StringC)vanilla + (*it);
		StringC outfile = (StringC)output + (*it);
		SampleC<SArray1dC<Point2dC> > pts = ReadPointsFromVectorFile(infile);
		//cout<<"Generating Sample of vectors"<<endl;
		SampleC<VectorC> pcapts;
		for(SampleIterC<SArray1dC<Point2dC> > it2(pts); it2; it2++)		
		{
			if(it2->Size() == 11)
				pcapts.Append(pca.PerformICAProjection(*it2).Copy());
		}
		cout<<"Writing "<<outfile<<endl;
		WriteVectorsToTSVFile(pcapts,outfile,infile);
	}

	return 0;
}

void WriteVectorsToTSVFile(SampleC<VectorC> const &pts, FilenameC const &file, FilenameC const &origfile)
{
	if(pts.Size() == 0)
	{
		//Just copy the original file to where it needs to go
		StringC cmd = "cp "+ (StringC) origfile+" "+(StringC) file;
		system(cmd);
	}
	else
	{
		UIntT vect_size = pts[0].Size();
		ofstream wfile;
		wfile.open(file.chars());
		wfile << pts.Size() << "\t" << vect_size<<endl;
		for(SampleIterC<VectorC> it(pts); it; it++)
		{
			for(UIntT i = 0; i < vect_size; i++)
				wfile<<(*it)[i]<<"\t";
			wfile<<"\n";
		}
		wfile.close();	
	}
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
		{
			cout<<"This file has vectors less than 22 in length: "<<file<<endl;
			// In this case no use ful thing will be returned
		}
		else
		{
			for(UIntT i = 0; i < numel; i++)
			{
				VectorC vect(vectsize); vect.Fill(0.0);
				for(UIntT j = 0; j < vectsize; j++)
				{
					myfile >> vect[j];
				}
				points.Append(VectorToPoints(vect).Copy());
			}
		}
		myfile.close();
	}
	return points;
}
