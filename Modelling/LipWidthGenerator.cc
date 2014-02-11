//ShapeModel Code
// Author - Bud Goswami
// Date - 9/11/09

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
#include "BSplineC.hh"
#include "ModellingUtilityFunctions.hh"
using namespace RavlN;
using namespace RavlImageN;

int main(int nargs,char **argv) 
{
	//CREATE THE REQUIRED INPUT OPTIONS
	OptionC opt(nargs,argv);
	//IO Options
	FilenameC lip = opt.String("l","in.txt","Text File with Lip-Co-ordinates");
	opt.Check();
	
	//Load the cascades
	DListC<SArray1dC<Point2dC> > lp_pts = LoadTRESData(lip);
	//This program will simply iterate through all the loaded tres files and compute the lip width
	for(DLIterC<SArray1dC<Point2dC> > it(lp_pts); it; it++)
	{
		cout<<Sqrt((*it)[8].SqrEuclidDistance((*it)[0]))<<endl;
	}
	return 0;	
}
