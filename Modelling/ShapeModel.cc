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

using namespace RavlN;
using namespace RavlImageN;

DListC<SArray1dC<Point2dC> > LoadTRESData(const FilenameC& f);
SArray1dC<Point2dC> ReorderPoints(const SArray1dC<Point2dC> &lp_pts); //Reorder Spline Points for EJ Landmarked Data
Array1dC<Point2dC> SArrayToArray(const SArray1dC<Point2dC> &pts);
ImageC<RealRGBValueC> LoadImage(const FilenameC& f);
VectorC ArrayToVector(const Array1dC<Point2dC> &arr);
Array1dC<Point2dC> PerformNormalisation(const Array1dC<Point2dC> &arr,const UIntT &n_type);
Array1dC<Point2dC> RotationNormalise(const Array1dC<Point2dC> &arr);
Array1dC<Point2dC> TranslationNormalise(const Array1dC<Point2dC> &arr);
Array1dC<Point2dC> PerformProjection(const Affine2dC &aff, const Array1dC<Point2dC> &pts);
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
	IntT ntype = opt.Int("t",0,"Normalisation type - 0-Raw data,1-translation,2-translation+rotation,3-translation+rotation+scale");
	opt.Check();
	
	//Load the cascades
	DListC<SArray1dC<Point2dC> > lp_pts = LoadTRESData(lip);
	//Load the images
	DListC<StringC> lp_imgs = img_dir.FiltList("*.ppm");
	//Now iterate through both the lists and draw the resulting lip points as well as the B-spline control points
	DLIterC<SArray1dC<Point2dC> > pts_it(lp_pts);
	//Marker Colours:
	RealRGBValueC tres(255,10,10),cps(10,255,10),spl(10,10,255);
	//cout<<"There are "<<lp_imgs.Size()<<" image frames and "<<lp_pts.Size()<<" data frames"<<endl;
	//IntT framenum = 0;
	for(DLIterC<StringC> img_it(lp_imgs); pts_it;img_it++)
	{
		//For the given SArray1dC<Point2dC> obtain the B-Spline Control Points
		BSplineC bspl(order,ncp, BSplineC::UOPEN);
		Array1dC<Point2dC> lip_cp = bspl.CalculateControlPoints(SArrayToArray((*pts_it)),BSplineC::UOPEN,BSplineC::CHORDLENGTH);
		//Now draw both of these into the image
		FilenameC im_name = img_dir + (*img_it);
		ImageC<RealRGBValueC> img_draw((LoadImage(im_name)).Copy());
		//cout<<"Rendering curve"<<endl;
		Array1dC<Point2dC> spl_pts = bspl.RenderCurve(resolution);
		pts_it++;
		//Perform Normalisation
		Array1dC<Point2dC> normed_cp = PerformNormalisation(lip_cp,ntype);
		VectorC cps_vec = ArrayToVector(normed_cp);
		//Output the control Points
		for(UIntT i = 0; i < cps_vec.Size(); i++)
		{
			cout<<cps_vec[i]<<"\t";
		}
		cout<<endl;
		/*for(Array1dIterC<Point2dC> it(lip_cp); it; it++)
		{
				//Draw B-Spline unnormalised control points on image
				DrawCross(img_draw,cps,(*it),4);		
		}
		for(SArray1dIterC<Point2dC> it((*pts_it)); it; it++)
		{
			//Draw the exrtacted tres points onto image
				DrawCross(img_draw,tres,(*it),2);
		}
		for(Array1dIterC<Point2dC> it(spl_pts); it; it++)
		{
			//Draw the rendered B-Spline curve on to image
			DrawCross(img_draw,spl,(*it),1);		
		}
		for(Array1dIterC<Point2dC> it(normed_cp); it; it++)
		{
			//Draw the normalised data curve on to image
			DrawCross(img_draw,spl,(*it),4);		
		}
		if(!Save("@X:Image Frame",img_draw)) cerr<<"Could not show the resulting spline overlaid image"<<endl;
		*/
	}
	return 0;	
}

DListC<SArray1dC<Point2dC> > LoadTRESData(const FilenameC& f)
//: Generate a data structure that contains the tres specified lip-shape co-ordinates
{
	//LOAD THE LANDMARKED DATA
	DListC<SArray1dC<Point2dC> > lp_coords;
	IStreamC is(f);
	while(!is.IsEndOfStream())
	{
		SArray1dC<Point2dC> pts;
		is >> pts;
		if(pts.Size() > 0)
		{
			//If data used is EJ, need to reverse order
			SArray1dC<Point2dC> to_append = ReorderPoints(pts);
			lp_coords.Append(to_append);			
		}				
	}	
	return lp_coords;
}

SArray1dC<Point2dC> ReorderPoints(const SArray1dC<Point2dC> &lp_pts)
//:This method is specifically applicable to the data provided by EJ using LP to landmark lip points
//:Since his method uses the CImg library, his x and y coordinates are reversed
//:EJ data is loaded straight into a Point2dC object, therefore, Old.Row() = New.Col() and vice versa.
{
	SArray1dC<Point2dC> new_lp(lp_pts.Size());
	for(UIntT i = 0; i < lp_pts.Size(); i++)
	{
		new_lp[i] = Point2dC(lp_pts[i].Col(), lp_pts[i].Row());
	}
	return new_lp;
}
 
Array1dC<Point2dC> SArrayToArray(const SArray1dC<Point2dC> &pts)
//:Data structure converter
{
	Array1dC<Point2dC> res(pts.Size());
	for(UIntT i = 0; i < pts.Size(); i++)
	{
		res[i] = pts[i].Copy();
	}
	return res;
}

ImageC<RealRGBValueC> LoadImage(const FilenameC& f)
//:Load Image
{
	ImageC<RealRGBValueC> img;
	if(!Load(f,img)) cerr<<"Could not load the image file "<<f<<endl;
	return img;
}

VectorC ArrayToVector(const Array1dC<Point2dC> &arr)
//:Convert the control point array to a vector - {X0....Xm,Y0...Ym} i.e. {Col[0]....Col[m],Row[0]...Row[m]}
{
	VectorC vec(arr.Size()*2);
	vec.Fill(0);
	for(UIntT i = 0 ; i < arr.Size() ; i++)
	{
		vec[i] = arr[i].Col();
		vec[i + arr.Size()] = arr[i].Row();
	}
	return vec;
}

Array1dC<Point2dC> TranslationNormalise(const Array1dC<Point2dC> &arr)
//:Perform translation normalisation
{
	//Compute the mean of the points
	Point2dC mean(0.0,0.0);
	mean = arr.Sum();
	mean = mean / (RealT)(arr.Size());
	//No subtract the mean position from each element inside the Array1dC
	Array1dC<Point2dC> res = arr - mean;
	return res;
}

Array1dC<Point2dC> RotationNormalise(const Array1dC<Point2dC> &arr)
//:	Perform rotation and translation normalisation using the lip corners
{
	//Lip Corners 
	Tuple2C<Point2dC,Point2dC> lip_corners(arr[0].Copy(),arr[8].Copy());
	//Horizontal x-axis
	Point2dC lpt(lip_corners.Data1().Copy());Point2dC rpt(lip_corners.Data1().Row(),lip_corners.Data2().Col());
	//Compute rotation angle
	Vector2dC x_axis = rpt - lpt;
	Vector2dC lip_crn = lip_corners.Data2() - lip_corners.Data1();
	RealT rot = lip_crn.Angle() - x_axis.Angle();    
	Vector2dC translation(0.0,0.0); //Translation is the mean vector
	RealT skew = 0.0;
	Vector2dC scaling(1.0,1.0);
	Affine2dC affine;
	affine.Compose(translation,scaling,skew,rot);
	//Apply the affine transformation to our points
	Array1dC<Point2dC> res = PerformProjection(affine,arr);
	return res;
}

Array1dC<Point2dC> PerformNormalisation(const Array1dC<Point2dC> &arr,const UIntT &n_type)
{
	Array1dC<Point2dC> norm_data;
	switch(n_type)
	{
		case 0:
		{
			//RAW DATA
			norm_data = arr.Copy();
			break;
		}
		case 1:
		{
			//TRANSLATION NORMALISE
			norm_data = TranslationNormalise(arr);
			break;
		}
		case 2:
		{
			//ROTATION NORMALISE
			Array1dC<Point2dC> rnorm_data = RotationNormalise(arr);
			norm_data = TranslationNormalise(arr);
			break;
		}
		case 3:
		{
			//SCALE NORMALISE
			Array1dC<Point2dC> rnorm_data = RotationNormalise(arr);
			norm_data = TranslationNormalise(arr);
		}
	}
	return norm_data;
}

Array1dC<Point2dC> PerformProjection(const Affine2dC &aff, const Array1dC<Point2dC> &pts)
{
	Array1dC<Point2dC> res(pts.Size());
	for(UIntT i = 0; i < pts.Size(); i++)
	{
		Vector2dC vec(pts[i].Row(),pts[i].Col());
		Vector2dC proj_vec = aff * vec;
		Point2dC proj_pt(proj_vec.Row(),proj_vec.Col());
		res[i] = proj_pt.Copy();
	}
	return res;
}
