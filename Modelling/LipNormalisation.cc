//GetBSplineControlPoints.cc
//Bud Goswami
//Purpose - given a tres file, use Mathematica B-Spline Interface to obtain the B-Spline Control Points for a heap of frames. Following this, save the results into a SampleC<SArray1dC<Point2dC> > output stream

//Required Libs
#include "Ravl/Affine2d.hh"
#include "Ravl/Array2d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/SArray2d.hh"
#include "Ravl/SArray1dIter.hh"
#include "Ravl/DArray1d.hh"
#include "Ravl/DP/SequenceIO.hh"
#include "Ravl/DP/FileFormatIO.hh"
#include "Ravl/IO.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/DrawCircle.hh"
#include "Ravl/Image/DrawLine.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/LinePP2d.hh"
#include "Ravl/Math.hh"
#include "Ravl/Matrix.hh"
#include "Ravl/Option.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/SArray1d.hh"
#include "Ravl/SArray2d.hh"
#include "Ravl/SArray1dIter.hh"
#include "Ravl/Stream.hh"
#include "Ravl/Tuple2.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Vector2d.hh"
#include "Ravl/Vector.hh"
#include "Ravl/Image/WarpAffine.hh"
#include <iostream>
#include <fstream>
//#include "BSplineC.hh"
#include "MBspline.hh"
#include "ModellingUtilityFunctions.hh"
#include "UtilityFunctions.hh"
#include "Align.hh"
#include <cstdio>
#include "ASMRotationNormalisation.hh"
#include "ASMAffineNormalisation.hh"
#include "Ravl/PatternRec/SampleIter.hh"
#include "Ravl/Image/AAMAppearance.hh"
using namespace RavlN;
using namespace RavlImageN;


 
int main(int nargs,char **argv) 
{
	//CREATE THE REQUIRED INPUT OPTIONS
	OptionC opt(nargs,argv);
	//IO Options
	FilenameC lip = opt.String("l","in.spl","Text File with Lip-Spline-Co-ordinates");
	opt.Check();
	
	ModellingUtilityFunctions mufc;
	//(1)Load the data
	SampleC<SArray1dC<Point2dC> > lip_feats;
	IStreamC myfile(lip);
	myfile>>lip_feats;
	//cout<<"Number of lip features is : "<<lip_feats.Size()<<endl;
	for(SampleIterC<SArray1dC<Point2dC> > it(lip_feats); it; it++)
	{
		//~ Point2dC center = ComputeCentroid(*it).Copy();
		//~ for(SArray1dIterC<Point2dC> it2(*it); it2; it2++)
		//~ {
			//~ (*it2) -= center;
		//~ }
		VectorC vect(PointsToVector(*it));
		for(UIntT i = 0; i < vect.Size(); i++)
			cout<<vect[i]<<"\t";
		cout<<endl;
	}
	return 0;	
}
