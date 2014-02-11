//File: LipInterpolation.cc
//Input: Query Image File, Type of Colour Space, Type of Cost Function
//Process: Uses Cascaded MCD to perform automatic lip segmentation and then fits a B-Spline Curve along the outer lip contour
//Author: Bud Goswami
//Date: 07.04.09
#include "Ravl/Array1d.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/Boundary.hh"
#include "Ravl/Crack.hh"
#include "Ravl/DP/SequenceIO.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Ellipse2d.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/DrawEllipse.hh"
#include "Ravl/Image/DrawFrame.hh"
#include "Ravl/Image/DrawCircle.hh"
#include "Ravl/Image/OpenCVConvert.hh"
#include "Ravl/Image/ImageRectangle.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/IO.hh"
#include "Ravl/MeanCovariance.hh"
#include "Ravl/Option.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/SArray1d.hh"
#include "Ravl/StdConst.hh"
#include "Ravl/String.hh"
#include "Ravl/SumsNd2.hh"
#include "opencv/highgui.h"
#include "opencv/cv.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <ctype.h>
#include "Omni/ColossusException.hh"
#include "Omni/FaceScan.hh"
#include "Omni/DetectFace.hh"
#include "Omni/FaceDetectionModel.hh"
#include "Omni/DetectedFace.hh"
#include "Omni/RawFace.hh"
#include "Omni/DisplayFace.hh"
//BUD SW
#include "ColourConvert.hh"//Contains colour conversion functions
#include "MouthRegion.hh"	//Contains mouth-region identification Functions 
#include "CostFunctions.hh"//Contains functions related to cluster grouping
#include "CascadedMCD.hh"
#include "BSplineC.hh"
using namespace RavlImageN;
using namespace RavlConstN;
using namespace OmniN;

template <typename T>
DListC<Index2dC> ClusterIdent(DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > &finalpop, const ImageC<T> &img, const IntT &segtype);
//Lip Boundary Search Functions
//Given a set of lip points, get best fitting ellipse
Ellipse2dC FitLipEllipse(const DListC<Index2dC> &lp);
//Obtain Points from the ellipse
Array1dC<Point2dC> GetEllipsePoints(const Ellipse2dC &ell, const UIntT &res);
//Given a list of lip pixels, obtain the mean covariance
template <typename T>
MeanCovarianceC ComputeMeanCovariance(const DListC<Index2dC> &arr, const ImageC<T> &img);
//Get Lip Bounary Points Using the J-Criterion Function
template <typename T>
Array1dC<Point2dC> GetLipBoundaryPointsJ(const BSplineC &bs, const UIntT &res, const ImageC<T> &img, const MeanCovarianceC &sk, const MeanCovarianceC &lp);
//Compute the J Functions
RealT ComputeCriterionFunctionJ(const VectorC &vec, const MeanCovarianceC &sk, const MeanCovarianceC &lp);
//Given a set of values, get the index of the last Zero Crossing
IntT GetLastZeroCrossing(const Array1dC<RealT> &pts);
//Get An Image Point From Any Image
template <typename T>
VectorC GetImagePoint(const ImageC<T> &img, const Point2dC &pt);
//Function to do all the lip boundary search stuff
template <typename T>
void DoLipBoundarySearch(const ImageC<T> &img, const ImageC<RealRGBValueC> &orig, const MeanCovarianceC &sk, const MeanCovarianceC &lp, const UIntT &ord, const UIntT &nc, const Array1dC<Index2dC> lip_pix, const ImageRectangleC &bfh, const UIntT &resolution);

int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	DirectoryC qimg_dir = opt.String("i","TestSet/","Input Query Image Directory");
	FilenameC cascade_name = opt.String("c","/usr/share/opencv/haarcascades/haarcascade_frontalface_alt.xml","Path for the Haar Classifier");
	UIntT pix_type = opt.Int("p",0,"Pixel Output Type - 0=RealRGBValueC, 1=NormalisedRG, 2=CIELab, 3 = PseudoHue");
	IntT segtype = opt.Int("t",0,"Type of Cluster Identification 0 - DensityCG, 1 - JCC");
	RealT h = opt.Real("h",0.75,"h-value for MCD algorithm"); 
	//B-Spline Parameters
	IntT order = opt.Int("d",3,"Spline Order");
	IntT ncp = opt.Int("n",11,"Number of Control Points Required");
	IntT resolution = opt.Int("r",10,"Number of Rendered Points Required");	
	opt.Check();
	
	DPIPortC<ImageC<RealRGBValueC>  > in;
   if(!OpenISequence(in,qimg)) { cerr<<"Could not load file "<<qimg<<endl;}
   ImageC<RealRGBValueC> src;
   //Also do the directoryC stuff!!
   while(in.Get(src)) 
   {
		//Get location of mouth-region ROI
		IndexC window = 50;
		cout<<"OBTAINING MOUTH REGION ROI"<<endl;
		ImageRectangleC botfacehalf = GetMouthRegion(qimg, cascade_name, window);
		DListC<Index2dC> lip_pix;
		//Now we have the necessary mouth-region window, we can perform colour conversion and clustering
  		switch(pix_type)
  		{
  			case 0:
  			{
  				ImageC<RealRGBValueC> m_img(GetRGBVals(qimg),botfacehalf);
  				if(!Save("@X:Mouth-region ROI",m_img)) cerr<<"Could not show mouth region roi"<<endl;
  				CascadedMCD mcd(m_img,h);
  				cout<<"PERFORMING LIP SEGMENTATION"<<endl;
  				DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > finalpop = mcd.Apply();
  				lip_pix = ClusterIdent(finalpop, m_img, segtype);
  				MeanCovarianceC lipmc = ComputeMeanCovariance(lip_pix, m_img);
  				DoLipBoundarySearch(m_img, src, finalpop.First().Data1(), lipmc,order,ncp,lip_pix, botfacehalf,resolution);
  				break;
  			}
 			case 1:
  			{
  				ImageC<VectorC> m_img(GetNormRGVals(qimg),botfacehalf);
  				CascadedMCD mcd(m_img,h);
   			cout<<"PERFORMING LIP SEGMENTATION"<<endl;
   			DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > finalpop = mcd.Apply();
  				lip_pix = ClusterIdent(finalpop, m_img, segtype);
 				MeanCovarianceC lipmc = ComputeMeanCovariance(lip_pix, m_img);
  				DoLipBoundarySearch(m_img, src, finalpop.First().Data1(), lipmc,order,ncp,lip_pix, botfacehalf,resolution);
  				break;
  			}
 			case 2:
  			{
  				ImageC<VectorC> m_img(GetCIELabVals(qimg),botfacehalf);
  				CascadedMCD mcd(m_img,h);
  				cout<<"PERFORMING LIP SEGMENTATION"<<endl;
  				DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > finalpop = mcd.Apply();
  				lip_pix = ClusterIdent(finalpop, m_img, segtype);
  				MeanCovarianceC lipmc = ComputeMeanCovariance(lip_pix, m_img);
  				DoLipBoundarySearch(m_img, src, finalpop.First().Data1(), lipmc,order,ncp,lip_pix, botfacehalf,resolution);
  				break;
  			}
 			case 3:
  			{
  				ImageC<VectorC> m_img(GetPseudoHueVals(qimg),botfacehalf);
  				CascadedMCD mcd(m_img,h);
  				cout<<"PERFORMING LIP SEGMENTATION"<<endl;
  				DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > finalpop = mcd.Apply();
  				lip_pix = ClusterIdent(finalpop, m_img, segtype);
  				MeanCovarianceC lipmc = ComputeMeanCovariance(lip_pix, m_img);
  				DoLipBoundarySearch(m_img, src, finalpop.First().Data1(), lipmc,order,ncp,lip_pix, botfacehalf,resolution);
  				break;
  			}
  		}	
   }
	return 0;
}

template <typename T>
DListC<Index2dC> ClusterIdent(DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > &finalpop, const ImageC<T> &img, const IntT &segtype)
{
	DListC<Index2dC> lip;
	//Method just chooses the type of function to use for cluster identification depending on segtype
	switch(segtype)
	{
		case 0:
		{
			lip = DensityCG(finalpop, img);
			break;
		}
		case 1:
		{
			lip = JCC(finalpop, img);
			break;
		}
	}
	return lip;
}

Ellipse2dC FitLipEllipse(const DListC<Index2dC> &lp)
{
	SArray1dC<Point2dC> lip_pts(lp.Size());
	UIntT i = 0;
	for(DLIterC<Index2dC> it(lp); it; it++)
		lip_pts[i++] = Index2dC((*it).Row(), (*it).Col());
	Ellipse2dC res;
	FitEllipse(lip_pts,res);
	return res;
}

Array1dC<Point2dC> GetEllipsePoints(const Ellipse2dC &ell, const UIntT &res)
{
	Array1dC<Point2dC> out(res);
	RealT inc = (2.0*pi)/(RealT)res;
	RealT angle = 0.0;
	for(UIntT i = 0; i < out.Size(); i++)
	{
		out[i] = ell.Point(angle);
		angle += inc;
	} 
	return out;
}

template <typename T>
MeanCovarianceC ComputeMeanCovariance(const DListC<Index2dC> &arr, const ImageC<T> &img)
{ 
	SizeT dim = img[img.Frame().TopLeft()].Size();
	SumsNd2C sum(dim);
	for(DLIterC<Index2dC> it(arr); it; it++)
	{
		VectorC v(dim);
		for(UIntT i = 0; i < dim; i++)
			v[i] = img[*it][i];
		sum += v.Copy();
	}
	return sum.MeanCovariance().Copy();
}

template <typename T>
Array1dC<Point2dC> GetLipBoundaryPointsJ(const BSplineC &bs, const UIntT &res, const ImageC<T> &img, const MeanCovarianceC &sk, const MeanCovarianceC &lp)
{
	Array1dC<Point2dC> out(res);
	UIntT arr_ind = 0;
	for(RealT pval= 0.0; pval < 1.0; pval+= (1.0/(RealT)res))		
	{
		LinePP2dC normal = bs.GetNormalLine(pval);
		Array1dC<RealT> scores(30);
		for(IntT i = -10; i <= 20; i++)
		{
			Point2dC cur_pt = normal.Point(i);
			VectorC vec = GetImagePoint(img, cur_pt);
			scores[i+10] = ComputeCriterionFunctionJ(vec, sk, lp);
		}
		//Now we have an array of data, compute get last zero crossing
		out[arr_ind] = scores[GetLastZeroCrossing(scores)];
		arr_ind++;
	}
	return out;
}

RealT ComputeCriterionFunctionJ(const VectorC &vec, const MeanCovarianceC &sk, const MeanCovarianceC &lp)
{
	RealT jskin = -Log(sk.Covariance().Det()) - sk.MahalanobisDistance(vec);
	RealT jlip = -Log(lp.Covariance().Det()) - lp.MahalanobisDistance(vec);
	return (jlip - jskin);
}
IntT GetLastZeroCrossing(const Array1dC<RealT> &pts)
{
	//Since we are going from -10 to +20
	IntT pt = 10;
	for(IntT i = pts.Size() - 1; i > 1; i--)
	{
	 	if((pts[i] < 0)&&( pts[i-1] > 0))
	 		pt = i - 1;
	}
	return pt;
}

template <typename T>
VectorC GetImagePoint(const ImageC<T> &img, const Point2dC &pt)
{
	UIntT dim = img[img.Frame().TopLeft()].Size();
	VectorC v(dim);
	Index2dC ind(pt.Row(), pt.Col());
	for(UIntT i = 0; i < dim; i++)
		v[i] = img[ind][i];
	return v;
}

template <typename T>
void DoLipBoundarySearch(const ImageC<T> &img, const ImageC<RealRGBValueC> &orig, const MeanCovarianceC &sk, const MeanCovarianceC &lp, const UIntT &ord, const UIntT &nc, const Array1dC<Index2dC> lip_pix, const ImageRectangleC &bfh, const UIntT &resolution)
{
	RealRGBValueC red(250,10,10); RealRGBValueC green(10,250,10); RealRGBValueC black(10,10,10);	RealRGBValueC blue(0,0,255);
	ImageC<UIntT> lip_bin = GetBinaryLipImage(lip_pix,bfh);
   if(!Save("@X:LB Image",lip_bin)) cerr<<"Could not save Lip Binary Image"<<endl;
 	//URS Bootstrap Lip Boundary Search Stuff
 	//GET ELLIPSE FROM LIP POINTS AND GET ELLIPSE CURVE POINTS
	Ellipse2dC ell=FitLipEllipse(lip_pix);
	bool EllipseFill = false;
	DrawEllipse(orig,blue,ell,EllipseFill);  
	if(!Save("@X:Ellipse Image",orig)) cerr<<"Could not save Lip Cracks Image"<<endl;
	//Fitting B-Spline to Ellipse
	BSplineC bspl(ord,nc, BSplineC::UOPEN);
	Array1dC<Point2dC> lp = GetEllipsePoints(ell, resolution);
	for(UIntT i = 0; i < lp.Size(); i++)
	{
	   DrawCircle(orig,red,lp[i],5);    	
	}
	Array1dC<Point2dC> cpts = bspl.CalculateControlPoints(lp, BSplineC::CHORDLENGTH, BSplineC::UOPEN);
	Array1dC<Point2dC> spl_pts = bspl.RenderCurve(resolution);
	for(Array1dIterC<Point2dC> it(cpts); it; it++)
	{
		DrawCross(orig,green,(*it),5);    
	}
	for(Array1dIterC<Point2dC> it(spl_pts); it; it++)
	{
		DrawCross(orig,blue,(*it),2);    
	}	
	if(!Save("@X: OutputImage ",orig)) cerr<<"Could Not Save Output Image"<<endl;   
	//Now perform boundary search stuff
	Array1dC<Point2dC> newlip = GetLipBoundaryPointsJ(bspl,resolution,img, sk,lp);
	//Given new boundary points, reparameterise the B-Spline
	cpts = bspl.CalculateControlPoints(newlip, BSplineC::CHORDLENGTH, BSplineC::UOPEN);
	spl_pts = bspl.RenderCurve(resolution);	
	for(Array1dIterC<Point2dC> it(spl_pts); it; it++)
	{
		DrawCross(src,black,(*it),2);    
	}	
}

