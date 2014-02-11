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
#include <cstdio>
using namespace RavlN;
using namespace RavlImageN;
Array1dC<Point2dC> PerformNormalisation(const Array1dC<Point2dC> &arr, const RealT &mode_lip, const UIntT &n_type);
void DrawArray(ImageC<RealRGBValueC> &img, const Array1dC<Point2dC> &pts, const RealRGBValueC &rgb);
int main(int nargs,char **argv) 
{
	//CREATE THE REQUIRED INPUT OPTIONS
	OptionC opt(nargs,argv);
	//IO Options
	FilenameC lip = opt.String("l","in.txt","Text File with Lip-Co-ordinates");
	DirectoryC img_dir = opt.String("i","in.avi","Input stream. ");//specify input file
	IntT order = opt.Int("d",3,"Spline Order");
	IntT ncp = opt.Int("n",11,"Number of Control Points Required");
	IntT resolution = opt.Int("r", 30, "Resolution for the rendered curve");
	IntT ntype = opt.Int("t",3,"Normalisation type - 0-Raw data,1-translation,2-translation+rotation,3-translation+rotation+scale");
	RealT mode_width = opt.Int("dl",82.1542,"Mode Lip Width");
	opt.Check();
	
	ModellingUtilityFunctions mufc;
	//(1)Load the data
	DListC<SArray1dC<Point2dC> > lp_pts = mufc.LoadTRESData(lip);
	//Load the images
	DListC<StringC> lp_imgs = img_dir.FiltList("*.ppm");
	//Now iterate through both the lists and draw the resulting lip points as well as the B-spline control points
	DLIterC<SArray1dC<Point2dC> > pts_it(lp_pts);
	//Marker Colours:
	RealRGBValueC tres(255,10,10),cps(10,255,10),spl(10,10,255);
	for(DLIterC<StringC> img_it(lp_imgs); pts_it;img_it++)
	{
		//(2)Given some data, perform normalisation throughout
		Array1dC<Point2dC> unnormalised_raw_data = mufc.SArrayToArray(*pts_it);
		Array1dC<Point2dC> normalised_raw_data = PerformNormalisation(unnormalised_raw_data, mode_width,ntype);
		//DEBUG STEP: Draw the points
		FilenameC im_name = img_dir + (*img_it);
		ImageC<RealRGBValueC> img_draw((mufc.LoadImage(im_name)).Copy());
		//DrawArray(img_draw,normalised_raw_data,tres);
		//DrawArray(img_draw,unnormalised_raw_data,cps);
		
		//(3)Perform parameterisation of this raw data using B-Splines
		BSplineC bspl(order,ncp, BSplineC::UOPEN);
		Array1dC<Point2dC> lip_cp = bspl.CalculateControlPoints(normalised_raw_data,BSplineC::UOPEN,BSplineC::CHORDLENGTH);	
		VectorC cps_vec = mufc.ArrayToVector(lip_cp);	
		for(UIntT i = 0; i < cps_vec.Size(); i++)
		{
			cout<<cps_vec[i]<<"\t";
		}
		cout<<endl;		
		pts_it++;//increase the points iterator
		//getchar();
	}
	return 0;	
}

Array1dC<Point2dC> PerformNormalisation(const Array1dC<Point2dC> &arr, const RealT &mode_lip, const UIntT &n_type)
{
	ModellingUtilityFunctions mu;
	RealT mean_width = 82.1542;
	Point2dC translation(0.0,0.0); Vector2dC scaling(1.0,1.0); RealT rot = 0.0;	
	Array1dC<Point2dC> normed_pts;
	switch(n_type)
	{
		case 0:
		{
			//RAW DATA
			break;
		}
		case 1:
		{
			//TRANSLATION NORMALISE
			translation = -mu.TranslationNormalise(arr);
			normed_pts =  mu.AffineNormalise(arr, translation, scaling,rot);
			break;
		}
		case 2:
		{
			//ROTATION NORMALISE
			//First translate normalise
			translation = -mu.TranslationNormalise(arr);
			Array1dC<Point2dC> dummy_normed_pts =  mu.AffineNormalise(arr, translation, scaling,rot);
			//Now apply rotation around the origin
			rot = -mu.RotationNormalise(dummy_normed_pts);//cout<<"Rotation angle = "<<rot<<endl;
			normed_pts =  mu.AffineNormalise(dummy_normed_pts, translation, scaling,0.0);
			break;
		}
		case 3:
		{
			//SCALE NORMALISE
			translation = -mu.TranslationNormalise(arr);
			Array1dC<Point2dC> dummy_normed_pts =  mu.AffineNormalise(arr, translation, scaling,rot);
			scaling = mu.ScaleNormalise(mean_width,mode_lip);
			rot = -mu.RotationNormalise(dummy_normed_pts); //cout<<"Rotation angle = "<<rot<<endl;			
			normed_pts =  mu.AffineNormalise(dummy_normed_pts, translation, scaling,0.0);			
			break;
		}
	}
	return normed_pts;
}
void DrawArray(ImageC<RealRGBValueC> &img, const Array1dC<Point2dC> &pts, const RealRGBValueC &rgb)
{
	for(Array1dIterC<Point2dC> it(pts); it; it++)
	{
		DrawCross(img,rgb,(*it),2);
	}
	if(!Save("@XA:Drawn Image",img)) cerr<<"Could not show the drawing"<<endl;
}
