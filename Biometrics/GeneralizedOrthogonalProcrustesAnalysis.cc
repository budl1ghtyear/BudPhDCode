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
#include	"Ravl/SArray1dIter2.hh"
#include 	"Ravl/PatternRec/Sample.hh"
#include 	"Ravl/PatternRec/SampleIter.hh"
#include 	"Ravl/Vector.hh"
#include 	"BANCAGT.hh"
#include 	"UtilityFunctions.hh"
#include 	"MBspline.hh"
using namespace RavlN;
using namespace std;
//MAIN PROGRAM
SampleC<SArray1dC<Point2dC> > InputFromRFile(FilenameC const &file);
Point2dC ComputeShapeCentroid(SArray1dC<Point2dC> const &pts);
SArray1dC<Point2dC> AlignShapeToOrigin(SArray1dC<Point2dC> const &pts);
SArray1dC<Point2dC> IsomorphicScaling(SArray1dC<Point2dC> const &pts);
Point2dC ComputeL2Norm(SArray1dC<Point2dC> const &pts);
int main(int nargs, char** nargv)
{	
	OptionC opt(nargs,nargv);
	FilenameC rawdata = opt.String("f","/vol/vssp/lip-tracking/StatisticalShapeModelsOfLips/TrainingData/C1Raw","Filenames for use in projection");
	FilenameC output = opt.String("o","/vol/vssp/lip-tracking/StatisticalShapeModelsOfLips/TrainingData/GPAC1Raw","Output");
	opt.Check();
	
	SampleC<SArray1dC<Point2dC> > data = InputFromRFile(rawdata);
	cout<<"For an element: "<<data[0]<<endl;
	cout<<"Centroid is: "<<ComputeShapeCentroid(data[0])<<endl;
	cout<<"Centroid normalisation is: "<<AlignShapeToOrigin(data[0])<<endl;
	cout<<"L2 Norm is : "<<ComputeL2Norm(AlignShapeToOrigin(data[0]))<<endl;
	cout<<"Isomorphic scaling is: \n"<<IsomorphicScaling(AlignShapeToOrigin(data[0]))<<endl;
	return 0;
}

SampleC<SArray1dC<Point2dC> > InputFromRFile(FilenameC const &file)
{
	SampleC<SArray1dC<Point2dC> > points;
	ifstream myfile;
	myfile.open(file,ifstream::in);
	if(myfile.is_open())
	{
		UIntT vectsize(22);
		while(!myfile.eof())
		{
			VectorC vect(vectsize); vect.Fill(0.0);
			for(UIntT j = 0; j < vectsize; j++)
			{
				myfile >> vect[j];
			}
			if(vect.SumOfAbs() != 0.)
				points.Append(VectorToPoints(vect).Copy());
		}
		myfile.close();
	}
	return points;
}

SArray1dC<Point2dC>  ComputeMeanShape(SampleC<SArray1dC<Point2dC> > const &arr)
{
	VectorC sum(22); sum.Fill(0.0);
	for(SampleIterC<SArray1dC<Point2dC> > it(arr); it ; it++ )
	{
		sum += PointsToVector(*it);
	}
	sum /= (RealT)(arr.Size());
	return VectorToPoints(sum);
}

Point2dC ComputeShapeCentroid(SArray1dC<Point2dC> const &pts)
{
	Point2dC sum(0.,0.);
	for(SArray1dIterC<Point2dC> it(pts); it; it++)
	{
		sum += (*it);
	} 
	return Point2dC(sum.Row()/(RealT)pts.Size(),sum.Col()/(RealT)pts.Size());
}

SArray1dC<Point2dC> AlignShapeToOrigin(SArray1dC<Point2dC> const &pts)
{
	Point2dC centroid = ComputeShapeCentroid(pts);
	SArray1dC<Point2dC> arr(pts.Size());
	for(SArray1dIter2C<Point2dC, Point2dC> it(pts,arr); it; it++)
	{
		it.Data2() = it.Data1() - centroid;
	}
	return arr;
}

Point2dC ComputeL2Norm(SArray1dC<Point2dC> const &pts)
{
	Point2dC sum(0.0,0.0);
	for(SArray1dIterC<Point2dC> it(pts); it; it++)
	{
		sum += Point2dC(Sqr(it->Row()),Sqr(it->Col())).Copy();
	}
	return Point2dC(Sqrt(sum.Row()),Sqrt(sum.Col()));
}

SArray1dC<Point2dC> IsomorphicScaling(SArray1dC<Point2dC> const &pts)
{
	Point2dC norm = ComputeL2Norm(pts);
	SArray1dC<Point2dC> arr(pts.Size());
	for(SArray1dIter2C<Point2dC,Point2dC> it(pts,arr); it; it++)
	{
		it.Data2() = Point2dC(it.Data1().Row()/norm.Row(),it.Data1().Col()/norm.Col());
	}
	return arr;
}
