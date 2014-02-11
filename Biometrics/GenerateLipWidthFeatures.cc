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
#include 	"Ravl/Ellipse2d.hh"
#include 	"Ravl/Pair.hh"
#include	"Ravl/LineABC2d.hh"
#include	"Ravl/SArray1dIter3.hh"
using namespace RavlN;
using namespace std;
//MAIN PROGRAM

void WriteVectorsToTSVFile(SampleC<VectorC> const &pts, FilenameC const &file, FilenameC const &origfile);
SampleC<SArray1dC<Point2dC> > ReadPointsFromVectorFile(FilenameC const &file);
DListC<StringC> GetFeatureFiles(DirectoryC const &dir, StringC const &ext, UIntT const &type=0 );
SampleC<VectorC> ComputeFrameDifferences(SampleC<SArray1dC<Point2dC> > const &pts);
SampleC<VectorC> ComputeCombinedFeatures(SampleC<SArray1dC<Point2dC> > const &shapes, SampleC<SArray1dC<Point2dC> > const &diffs);
PairC<Point2dC> GetLipCornerPoints(SArray1dC<Point2dC> const &lip_bound_pts);
int main(int nargs, char** nargv)
{	
	OptionC opt(nargs,nargv);
	DirectoryC shape_features = opt.String("s","/vol/vssp/lip-tracking/BiometricLipFiles/VanillaShapeFeatures/ShapeFeatures","Directory containing vanilla features");
	StringC ext = opt.String("e",".TSV","Extension of the vanilla features");
	FilenameC outwidthfile = opt.String("ow","/vol/vssp/lip-tracking/lipwidths.txt","Output file containing only lip widths");
	FilenameC outheightfile = opt.String("oh","/vol/vssp/lip-tracking/lipheights.txt","Output file containing only lip widths");
	UIntT feattype = opt.Int("t",0,"Type of feature - client (0) or world (1)");
	opt.Check();
	
	DListC<StringC> shape_feat_files = GetFeatureFiles(shape_features, ext, feattype);
	//~ cout<<"The file we are going to test is: "<<shape_feat_files.Nth(0);
	//~ SampleC<SArray1dC<Point2dC> > lipshapes = ReadPointsFromVectorFile(shape_feat_files.Nth(0));
	//~ PairC<Point2dC> extremes = GetLipCornerPoints(lipshapes.First());
	//~ cout<<"Lip Points are: \n"<<lipshapes.First()<<"\n Extremes: \n"<<extremes<<endl;
	SampleC<RealT> widths;
	SampleC<RealT> heights;
	RealT maxw(0.),maxh(0.);
	ofstream myoutfile; ofstream myhoutfile;
	myoutfile.open(outwidthfile.chars());myhoutfile.open(outheightfile.chars());
	FilenameC maxwfile,maxhfile;
	//~ //Now we have the list of files we need to input, we need to read them and perform our frame differencing operation on the file:
	for(UIntT i = 0; i < shape_feat_files.Size(); i++)
	{
		FilenameC orgshapefile = shape_feat_files.Nth(i);
		SampleC<SArray1dC<Point2dC> > lipshapes = ReadPointsFromVectorFile(orgshapefile);
		for(SampleIterC<SArray1dC<Point2dC> > it(lipshapes); it; it++)
		{
				PairC<Point2dC> extremes = GetLipCornerPoints(*it);
				RealT width = Abs(extremes.Data1().X() - extremes.Data2().X()); //Insert the current lip width into the sample collection
				RealT height = Abs(extremes.Data1().Y() - extremes.Data2().Y()); //Insert the current lip width into the sample collection
				myoutfile << width <<"\n";
				myhoutfile << height <<"\n";
				widths.Append(width);
				heights.Append(height);				
				if(width > maxw)
				{
					maxw = width;
					maxwfile = orgshapefile;
				}
				if(height > maxh)
				{
					maxh = height;
					maxhfile = orgshapefile;
				}
		}
	}
	myoutfile.close();
	myhoutfile.close();
	cout<<"Maximum width = "<<maxw<<"\t Maximum heigth = "<<maxh<<endl;
	cout<<" Width file: "<<maxwfile<<"\n Max Height File: "<<maxhfile<<endl;
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

PairC<Point2dC> GetLipCornerPoints(SArray1dC<Point2dC> const &lip_bound_pts)
{
	Point2dC centre;
	RealT major;
	RealT minor;
	RealT angle;
	Ellipse2dC lip_ellipse;
	FitEllipse(lip_bound_pts,lip_ellipse);
	lip_ellipse.EllipseParameters(centre, major, minor, angle);
	SArray1dC<RealT> major_axis_dist(lip_bound_pts.Size());
	major_axis_dist.Fill((RealT)100.0);
	SArray1dC<RealT> minor_axis_dist(lip_bound_pts.Size());
	minor_axis_dist.Fill((RealT)0.0);
	LineABC2dC major_line(-Sin(angle), Cos(angle), -(-Sin(angle)*centre.Row()+Cos(angle)*centre.Col()));
	LineABC2dC minor_line(Cos(angle), Sin(angle), -(Cos(angle)*centre.Row()+Sin(angle)*centre.Col()));
	for (SArray1dIter3C<Point2dC, RealT, RealT> it(lip_bound_pts, major_axis_dist, minor_axis_dist);it;it++)
	{
		it.Data2()=major_line.SqrEuclidDistance(it.Data1());
		it.Data3()=minor_line.SignedDistance(it.Data1());
	}
	IndexC leftBound = minor_axis_dist.IndexOfMin();
	if (major_axis_dist[leftBound] < 100)
	{
	}
	IndexC rightBound = minor_axis_dist.IndexOfMax();
	if (major_axis_dist[rightBound] < 100)
	{
	}
	return (PairC<Point2dC> (lip_bound_pts[leftBound].Copy(),lip_bound_pts[rightBound].Copy()));	
}
