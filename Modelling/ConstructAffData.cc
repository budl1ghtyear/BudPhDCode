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


void DrawArray(ImageC<RealRGBValueC> &img, const Array1dC<Point2dC> &pts, const RealRGBValueC &rgb);
int main(int nargs,char **argv) 
{
	//CREATE THE REQUIRED INPUT OPTIONS
	OptionC opt(nargs,argv);
	//IO Options
	FilenameC lip = opt.String("l","in.txt","Text File with Lip-Co-ordinates");
	//IntT order = opt.Int("d",3,"Spline Order");
	//IntT ncp = opt.Int("n",11,"Number of Control Points Required");
	RealT var = opt.Real("v",0.95,"Variation to be preserved in the parameter");
	//IntT resolution = opt.Int("r", 30, "Resolution for the rendered curve");
	IntT type = opt.Int("t",0,"Type of normalisation to be performed:0 - Perform 3 parameter GPA, 1 - Perform 4 parameter GPA, 2 - Perform 6 parameter GPA");
	opt.Check();
	
	ModellingUtilityFunctions mufc;
	//(1)Load the data
	DListC<SArray1dC<Point2dC> > lp_points = mufc.LoadTRESData(lip);
	DListC<SArray1dC<Point2dC> > lp_pts;
	SampleC<AAMAppearanceC> lip_feats;
	//UIntT counter = 0;
	for( DLIterC<SArray1dC<Point2dC> > pts_it(lp_points); pts_it; pts_it++)
	{
		ModellingUtilityFunctions mufc;
		//BSplineC bspl(order,ncp, BSplineC::UOPEN);
		//Array1dC<Point2dC> lip_cp = bspl.CalculateControlPoints(mufc.SArrayToArray((*pts_it)),BSplineC::UOPEN,BSplineC::CHORDLENGTH);
		//cout<<counter++<<endl;
		MBSplineC bspl;
		Array1dC<Point2dC> lip_cp = bspl.FitQuadraticSpline(mufc.SArrayToArray((*pts_it)));	
		lp_pts.Append(mufc.ArrayToSArray(lip_cp).Copy());
		lip_feats.Append(mufc.ArrayToSArray(lip_cp).Copy());
	}	
	
	//Now we have a list of points in a SampleC<AppearanceC>
	//AAMAppearanceMirrorC mirror;
	//mirror.Reflect(lip_feats);
	SampleVectorC smpl;
	switch(type)
	{
		case 0:
		{
			ASMNormalisationC asmdl(true);
			asmdl.Design(lip_feats, var,22);
			smpl = asmdl.GetNormalisedPoints();
			break;
		}
		case 1:
		{	
			ASMScaleRotationNormalisationC asmdl(true);
			asmdl.Design(lip_feats, var,22);
			smpl = asmdl.GetNormalisedPoints();
			break;
		}
		case 2:
		{
			ASMAffineNormalisationC asmdl(true);
			asmdl.Design(lip_feats, var,22);
			smpl = asmdl.GetNormalisedPoints();
			break;
		}
	}
	for(SampleIterC<VectorC> it(smpl); it; it++)
	{
		UIntT size = it->Size();
		for(UIntT i = 0; i < size; i++)
		{
			cout<<(*it)[i]<<"\t";
		}
		cout<<endl;
	}

	/*	
	DListC<SArray1dC<Point2dC> > gpa_pts;		
	switch(type)
	{	
		case 1:
		{
			//Perform Similarity Transform
			AlignC align;
			//~ DListC<SArray1dC<Point2dC> > dummy = align.PerformSimilarityAlignment(lp_pts,lp_pts.First());
			//~ gpa_pts=align.ComputeDeviation(dummy,lp_pts.First()).Copy();			
			gpa_pts = align.PerformSimilarityAlignment(lp_pts,lp_pts.First());
			break;
		}
		case 2:
		{
			AlignC align;
			//~ DListC<SArray1dC<Point2dC> > dummy = align.PerformAffineAlignment(lp_pts,lp_pts.First());
			//~ gpa_pts=align.ComputeDeviation(dummy,lp_pts.First()).Copy();
			gpa_pts = align.PerformAffineAlignment(lp_pts,lp_pts.First());
			break;
		}
	}	
	//Now iterate through both the lists and draw the resulting lip points as well as the B-spline control points
	for(DLIterC<SArray1dC<Point2dC> > it(gpa_pts); it; it++)
	{
		for(SArray1dIterC<Point2dC> sit(*it); sit; sit++)
		{
			cout<<(*sit).Col()<<"\t"<<(*sit).Row()<<"\t";
		}
		cout<<endl;
	}
	*/
	return 0;	
}
