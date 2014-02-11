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
	FilenameC lip = opt.String("l","in.txt","Text File with Lip-Co-ordinates");
	opt.Check();
	
	ModellingUtilityFunctions mufc;
	//(1)Load the data
	DListC<SArray1dC<Point2dC> > lp_points = mufc.LoadTRESData(lip);
	SampleC<SArray1dC<Point2dC> > lip_feats;
	//UIntT counter = 0;
	for( DLIterC<SArray1dC<Point2dC> > pts_it(lp_points); pts_it; pts_it++)
	{
		ModellingUtilityFunctions mufc;
		MBSplineC bspl;
		Array1dC<Point2dC> lip_cp = bspl.FitQuadraticSpline(mufc.SArrayToArray((*pts_it)));	
		lip_feats.Append(mufc.ArrayToSArray(lip_cp).Copy());
	}	
	//Now save the SampleC into an output stream
	FilenameC lip_file = lip.PathComponent()+"/" + lip.BaseNameComponent() + ".spl";
	OStreamC myfile(lip_file);
	myfile<<lip_feats;
	return 0;	
}
