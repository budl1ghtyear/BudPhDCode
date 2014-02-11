//LipClusteringC test program

//Just to test the Lip Clustering Test Class
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/SArray1d.hh"
#include "Ravl/SArray1dIter3.hh"
#include "Ravl/Option.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/String.hh"
#include "LipSegmentation.hh"
#include "ColourSpaceImage.hh"
#include "XMLToGTConv.hh"
#include "BANCAGT.hh"
#include "Ravl/Boundary.hh"
#include "Ravl/Image/MorphClose.hh"
#include "Ravl/Image/Erode.hh"
#include "Ravl/Ellipse2d.hh"
#include "Ravl/Image/DrawEllipse.hh"
#include "Ravl/Array2dIter3.hh"
#include "Ravl/OS/Date.hh"
#include "LipExtractorC.hh"
//#include "BSplineSolver.hh"
#include "BSplineC.hh"

//#define NUM_POINTS_EACH_SIDE 6
using namespace RavlN;
using namespace RavlConstN;
using namespace RavlImageN;

RealT ComputeSegmentationQuality(const ImageC<UIntT> &im, const ImageC<UIntT> &gt);
//Method to compute the image quality based on binary ground truth and segmented lip image

int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	FilenameC img_name = opt.String("i","Input.ppm","example input image");
	IntT cluster_type = opt.Int("ct",0,"ClusteringType - values allowed - CMCD=0,KSMCD=1,KM=2,FCM=3, default = CMCD");
	IntT pixel_type = opt.Int("pt",0,"Pixel Type - RGB=0, NormRG=1, HSV=2, YUV=3, CIE=4, PHue=5 , default = RGB=0");
	IntT label_type = opt.Int("lt",0,"Pixel labelling type - JNormal=0, JSkin=1, default =  JNormal=0 ");
	IntT region_type = opt.Int("rt",0,"Region Identification type - CC=0,CCentral=1, default =CC=0 ");
	RealT hvalue = opt.Real("h",0.75,"Confidence Value of MCD - values between 0.5 and 1.0, default = 0.75");
	UIntT numclust = opt.Int("n",3,"Number of clusters for KMeans or Fuzzy C means, default = 3, must be greater than 1");
	opt.Compulsory("i");
	opt.Check();
	
	//Load input image
	ImageC<RealRGBValueC> img;
	if(1==(!Load(img_name, img))) cerr<<"Could not load input image"<<endl;
	//if(!Save("@X:Original Image",img)) cerr<<"Could not show input image"<<endl;
	//cout<<"Instantiating class"<<endl;
	LipSegmentationC lipseg(img_name);
	Tuple3C<ImageC<UIntT>, MeanCovarianceC, MeanCovarianceC> result = lipseg.Apply((PixelType)pixel_type,(ClusteringType)cluster_type,(LabellingType) label_type,(RegionIdentificationType) region_type,hvalue, numclust);
	//if(!Save("@X:Binary Image",result.Data1())) cerr<<"Could not show binary lip file"<<endl; //Now we have the lip binary image in real space
	//Make full size Binary Lip Image
	ImageC<UIntT> lip_img(img.Frame(),0);
	DListC<Index2dC> lip_ind;
	for(Array2dIterC<UIntT> it(result.Data1()); it; it++)
	{
		Index2dC ind = it.Index();
		if((*it) == 255)
		{
			lip_img[ind] = 255;
			lip_ind.Append(ind.Copy());
		}
	}
	if(!Save("@X:Binary Lip Image",lip_img)) cerr<<"Could not save binary lip image"<<endl;

	LipExtractorC lext(lip_img,result.Data1().Frame(),20);
	Array1dC<Point2dC> lip_points(lext.GetLipBoundaryPointsArray());
	cout<<"Lip Points "<<lip_points<<endl; 
	//BSplineSolver bspl(lip_points, 4, 11, BSplineSolver::AVERAGING, BSplineSolver::CHORDLENGTH);
	BSplineC bspl(4, 11, BSplineC::UOPEN, false);
	Array1dC<Point2dC> spl_cp = bspl.CalculateControlPoints(lip_points,BSplineC::UOPEN,BSplineC::UNIFORM);
	//BSplineC bspl(4, 11, BSplineC::UPCLOSED, false);
	//Array1dC<Point2dC> spl_cp = bspl.CalculateControlPoints(lip_points,BSplineC::UPCLOSED,BSplineC::UNIFORM);
	cout<< "B-Spline Control Points Are - \n"<< spl_cp << endl;
	//Array1dC<Point2dC> spl_pts = bspl.RenderCurve(500);
	Array1dC<Point2dC> spl_pts = bspl.RenderCurve(100);
	ImageC<RealRGBValueC> spl_img(img.Copy());
	RealRGBValueC cp_col(20,255,20);RealRGBValueC pt_col(255,20,20); RealRGBValueC lip_col(20,20,20);
	
	for(Array1dIterC<Point2dC> it(spl_pts); it; it++)
	{
		DrawCross(spl_img,pt_col,Index2dC((*it).Row(),(*it).Col()),1);
	}
	
	for(Array1dIterC<Point2dC> it(spl_cp); it; it++)
	{
		DrawCross(spl_img,cp_col,Index2dC((*it).Row(),(*it).Col()),3);
	}
	for(Array1dIterC<Point2dC> it(lip_points); it; it++)
	{
		DrawCross(spl_img,lip_col,Index2dC((*it).Row(),(*it).Col()),4);
	}
	if(!Save("@X:Spline Lip Image",spl_img)) cerr<<"Could not save binary lip image"<<endl;
	/*
	BSplineC newbspl(4, bspl.GetKnotVector(),lip_points);
	cout<<"Points inserted into new B-Spline are - "<<lip_points<<endl;
	sArray1dC<Point2dC> newspl_pts = newbspl.RenderClosedCurve(100);
	ImageC<RealRGBValueC> newspl_img(img);
	for(Array1dIterC<Point2dC> it(newspl_pts); it; it++)
	{
		DrawCross(newspl_img,pt_col,Index2dC((*it).Row(),(*it).Col()),1);		
		//if(!Save("@X:New Spline Lip Image",newspl_img.Copy())) cerr<<"Could not save binary lip image"<<endl;
		//Sleep(0.1);
	}
	if(!Save("@X:New Spline Lip Image",newspl_img.Copy())) cerr<<"Could not save binary lip image"<<endl;
	
	cout<<"Newspl_pts size "<<newspl_pts.Size()<<endl;
	*/

	//BSPLINE STUFF
	/*BSplineC spl(3, 11,BSplineC::UPERIODIC);
	Array1dC<Point2dC> spl_cp = spl.CalculateControlPoints(lip_points, BSplineC::CHORDLENGTH, BSplineC::UOPEN, 3,11);
	cout<< "B Spline Control Points Extracted Are \n" << spl_cp << endl;
	Array1dC<Point2dC> spl_pts = spl.RenderCurve(30);
	cout<< "B Spline Curve Points Extracted Are \n" << spl_pts << endl;
	ImageC<RealRGBValueC> spl_img(img);
	RealRGBValueC cp_col(20,255,20);RealRGBValueC pt_col(255,20,20);
	
	for(Array1dIterC<Point2dC> it(spl_pts); it; it++)
	{
		DrawCross(spl_img,pt_col,Index2dC((*it).Row(),(*it).Col()),2);
	}
	
	for(Array1dIterC<Point2dC> it(spl_cp); it; it++)
	{
		DrawCross(spl_img,cp_col,Index2dC((*it).Row(),(*it).Col()),2);
	}
	if(!Save("@X:Spline Lip Image",spl_img)) cerr<<"Could not save binary lip image"<<endl;*/
	return 0;
}
