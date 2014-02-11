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
DListC<StringC> GetFeatureFiles(DirectoryC const &dir, StringC const &ext, UIntT const &type=0 );
SampleC<VectorC> ComputeFrameDifferences(SampleC<SArray1dC<Point2dC> > const &pts);
SampleC<VectorC> ComputeCombinedFeatures(SampleC<SArray1dC<Point2dC> > const &shapes, SampleC<SArray1dC<Point2dC> > const &diffs);
int main(int nargs, char** nargv)
{	
	OptionC opt(nargs,nargv);
	DirectoryC shape_features = opt.String("s","/vol/vssp/lip-tracking/BiometricLipFiles/VanillaShapeFeatures/ShapeFeatures","Directory containing vanilla features");
	DirectoryC diff_features = opt.String("d","/vol/vssp/lip-tracking/BiometricLipFiles/VanillaDiffFeatures/Features","Directory containing vanilla features");
	StringC ext = opt.String("e",".TSV","Extension of the vanilla features");
	DirectoryC output = opt.String("o","/vol/vssp/lip-tracking/BiometricLipFiles/VanillaDiffFeatures/Features/","Output directory in which to store the file");
	UIntT feattype = opt.Int("t",0,"Type of feature - client (0) or world (1)");
	opt.Check();
	
	DListC<StringC> shape_feat_files = GetFeatureFiles(shape_features, ext, feattype);
	DListC<StringC> diff_feat_files = GetFeatureFiles(diff_features, ext, feattype);
	
	//Now we have the list of files we need to input, we need to read them and perform our frame differencing operation on the file:
	for(UIntT i = 0; i < shape_feat_files.Size(); i++)
	{
		FilenameC orgshapefile = shape_feat_files.Nth(i);
		FilenameC orgdifffile;
		for(DLIterC<StringC> it(diff_feat_files); it; it++)
		{
			orgdifffile = (*it);
			if(orgdifffile.BaseNameComponent() == orgshapefile.BaseNameComponent())
				break;
		}
		FilenameC outfile = (StringC)output + "/" + (StringC)(orgshapefile.NameComponent());
		
		//cout<<"About to combine "<<orgshapefile<<"  With "<<orgdifffile<<endl;
		SampleC<VectorC> comb_vects = ComputeCombinedFeatures(ReadPointsFromVectorFile(orgshapefile),ReadPointsFromVectorFile(orgdifffile));
		WriteVectorsToTSVFile(comb_vects,outfile,orgshapefile);
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
		if(vectsize == 1)
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
DListC<StringC> GetFeatureFiles(DirectoryC const &dir, StringC const &ext, UIntT const &type)
//This method will return a list of feature files from within a directory depending on whether it is feature directory or world directory
{
	StringC extension = "*"+ext;
	
	DListC<StringC> files;
	if(type == 0)
	{
		//If normal features, files will be as: dir/<client_name>/*.TSV
		//Get all client dirs
		DListC<StringC> clients = dir.FiltList("*");
		//cout<<"Clients are : "<<clients<<endl;
		for(DLIterC<StringC> it(clients); it; it++)
		{
			DirectoryC curdir = (StringC)dir + "/" + (*it);
			DListC<StringC> fls = curdir.FiltList(extension);
			for(DLIterC<StringC> it2(fls); it2; it2++)
			{
				StringC filename = (StringC)dir + "/" + (*it) + "/" + (*it2);
				files.Append(filename);
			}
		} 
	}
	else if (type ==1)
	{
		//If world features, files will be as: dir/*.TSV
		DListC<StringC> fls = dir.FiltList(extension);
		for(DLIterC<StringC> it2(fls); it2; it2++)
		{
			StringC filename = (StringC)dir + "/" + (*it2);
			files.Append(filename);
		}		
	}
	else
	{
		cerr<<"SORRY INCORRECT FEATURE TYPE SPECIFIER. PLEASE ENTER 0 FOR NORMAL FEATURES and 1 FOR WORLD FEATURES"<<endl;
		exit(1);
	}
	return files;
}

SampleC<VectorC> ComputeFrameDifferences(SampleC<SArray1dC<Point2dC> > const &pts)
{
	SampleC<VectorC> res;
	if(pts.Size() == 0)
	{
		return res;
	}
	else
	{
		for(UIntT i = 1; i < pts.Size(); i++)
		{
			VectorC cur = PointsToVector(pts[i]).Copy();
			VectorC prev = PointsToVector(pts[i-1]).Copy();
			res.Append(cur - prev);
		}
		return res;
	}
}

SampleC<VectorC> ComputeCombinedFeatures(SampleC<SArray1dC<Point2dC> > const &shapes, SampleC<SArray1dC<Point2dC> > const &diffs)
{
	SampleC<VectorC> res;
	if(shapes.Size() == 0)
	{
		return res;
	}
	else
	{
		//Now simply combine the two vectors
		for(UIntT i = 1; i < shapes.Size(); i++)
		{
			VectorC shapevect = PointsToVector(shapes[i]).Copy();
			VectorC diffvect = PointsToVector(diffs[i-1]).Copy();
			res.Append(shapevect.Append(diffvect).Copy());
		}
		return res;
	}
}
