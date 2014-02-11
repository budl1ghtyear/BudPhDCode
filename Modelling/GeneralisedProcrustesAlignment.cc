//ConstructShapeModel Code
// Author - Bud Goswami
// Date - 9/11/09
// Process:
//1) Load the TRES DATA
//2) Perform either translation, rotation or scale normalisation on the TRES data
//3) Perform parameterisation of this raw data
//4) Output the data 

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
#include "ProcrustesAlign.hh"
#include "Align.hh"
#include <cstdio>
using namespace RavlN;
using namespace RavlImageN;


void DrawArray(ImageC<RealRGBValueC> &img, const Array1dC<Point2dC> &pts, const RealRGBValueC &rgb);
int main(int nargs,char **argv) 
{
	//CREATE THE REQUIRED INPUT OPTIONS
	OptionC opt(nargs,argv);
	//IO Options
	FilenameC lip = opt.String("l","in.txt","Text File with Lip-Co-ordinates");
	IntT order = opt.Int("d",3,"Spline Order");
	IntT ncp = opt.Int("n",11,"Number of Control Points Required");
	IntT resolution = opt.Int("r", 30, "Resolution for the rendered curve");
	opt.Check();
	
	ModellingUtilityFunctions mufc;
	//(1)Load the data
	DListC<SArray1dC<Point2dC> > lp_pts = mufc.LoadTRESData(lip);
	ProcrustesAlignC gpa(lp_pts.First());
	DListC<SArray1dC<Point2dC> > gpa_pts = gpa.Apply(lp_pts);

	//Now iterate through both the lists and draw the resulting lip points as well as the B-spline control points
	for( DLIterC<SArray1dC<Point2dC> > pts_it(lp_pts); pts_it; pts_it++)
	{
		ModellingUtilityFunctions mufc;
		BSplineC bspl(order,ncp, BSplineC::UOPEN);
		Array1dC<Point2dC> lip_cp = bspl.CalculateControlPoints(mufc.SArrayToArray((*pts_it)),BSplineC::UOPEN,BSplineC::CHORDLENGTH);
		VectorC cps_vec = mufc.ArrayToVector(lip_cp);
		//Output the control Points
		for(UIntT i = 0; i < cps_vec.Size(); i++)
		{
			cout<<cps_vec[i]<<"\t";
		}
		cout<<endl;	
	}
	return 0;	
}
