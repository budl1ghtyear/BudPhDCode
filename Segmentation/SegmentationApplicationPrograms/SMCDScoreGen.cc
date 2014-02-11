//File: BasicThesisSeg.cc
//Input: Query Image Directory, Ground Truth Directory, Type of Colour Space, Type of Cost Function
//Author: Bud Goswami
//Date: 18.03.09
#include "Ravl/Array1d.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/DrawFrame.hh"
#include "Ravl/Image/OpenCVConvert.hh"
#include "Ravl/Image/ImageRectangle.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/IO.hh"
#include "Ravl/Option.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/String.hh"
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
#include "MCD.hh"
using namespace RavlImageN;
using namespace OmniN;

template <typename T>
DListC<Index2dC> ClusterIdent(Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > &finalpop, const ImageC<T> &img, const IntT &segtype);

int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	DirectoryC qimg_dir = opt.String("i","TestSet/","Input Query Image Directory");
	DirectoryC gtqimg_dir = opt.String("g","GT/","Input Ground Truth Image Directory");
	FilenameC cascade_name = opt.String("c","/usr/share/opencv/haarcascades/haarcascade_frontalface_alt.xml","Path for the Haar Classifier");
	UIntT pix_type = opt.Int("p",0,"Pixel Output Type - 0=RealRGBValueC, 1=NormalisedRG, 2=RealHSVValueC, 3=RealYUVValueC, 4=CIELab, 5 = PseudoHue");
	IntT segtype = opt.Int("t",0,"Type of Cluster Identification 0 - DensityCG, 1 - JCC");
	RealT h = opt.Real("h",0.75,"h-value for MCD algorithm"); 
	opt.Check();
	
	DListC<StringC> qfile = qimg_dir.FiltList("*.ppm");
	for(DLIterC<StringC> it(qfile); it; it++)
	{
		FilenameC qimg = qimg_dir+(*it);
		FilenameC gtimg = gtqimg_dir + qimg.BaseNameComponent() + ".pbm";
	
		ImageC<RealRGBValueC> src;   
		if(!Load(qimg, src)) cerr<<"Loading RAVL Image Failed"<<endl;
		//Get location of mouth-region ROI
		IndexC window = 50;
		ImageRectangleC botfacehalf = GetMouthRegion(qimg, cascade_name, window);
		DListC<Index2dC> lip_pix;
		//Now we have the necessary mouth-region window, we can perform colour conversion and clustering
  		switch(pix_type)
  		{
  			case 0:
  			{
  				ImageC<RealRGBValueC> m_img(GetRGBVals(qimg),botfacehalf);
  				if(!Save("@X:Mouth-region ROI",m_img)) cerr<<"Could not show mouth region roi"<<endl;
  				MCD mcd(m_img,h);
  				Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > finalpop = mcd.Apply();
  				DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > lp;
  				lp.Append(finalpop);
  				lip_pix = ClusterIdent(lp, m_img, segtype);
  				break;
  			}
 			case 1:
  			{
  				ImageC<VectorC> m_img(GetNormRGVals(qimg),botfacehalf);
  				MCD mcd(m_img,h);
  				Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > finalpop = mcd.Apply();
  				DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > lp;
  				lp.Append(finalpop);
  				lip_pix = ClusterIdent(lp, m_img, segtype);
  				break;
  			}
  			case 2:
  			{
  				ImageC<RealHSVValueC> m_img(GetHSVVals(qimg),botfacehalf);
  				if(!Save("@X:Mouth-region ROI",m_img)) cerr<<"Could not show mouth region roi"<<endl;
  				MCD mcd(m_img,h);
  				Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > finalpop = mcd.Apply();
				lip_pix = ClusterIdent(finalpop, m_img, segtype);
   			break;
  			}
  			case 3:
  			{
  				ImageC<RealYUVValueC> m_img(GetYUVVals(qimg),botfacehalf);
  				if(!Save("@X:Mouth-region ROI",m_img)) cerr<<"Could not show mouth region roi"<<endl;
  				MCD mcd(m_img,h);
  				Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > finalpop = mcd.Apply();
				lip_pix = ClusterIdent(finalpop, m_img, segtype);
   			break;
  			}
 			case 4:
  			{
  				ImageC<VectorC> m_img(GetCIELabVals(qimg),botfacehalf);
  				MCD mcd(m_img,h);
  				Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > finalpop = mcd.Apply();
  				DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > lp;
  				lp.Append(finalpop);
  				lip_pix = ClusterIdent(lp, m_img, segtype);
  				break;
  			}
 			case 5:
  			{
  				ImageC<VectorC> m_img(GetPseudoHueVals(qimg),botfacehalf);
  				MCD mcd(m_img,h);
  				Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > finalpop = mcd.Apply();
  				DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > lp;
  				lp.Append(finalpop);
  				lip_pix = ClusterIdent(lp, m_img, segtype);
  				break;
  			}
  		}	
  		//Now we have the lip pixels as obtained by algorithm
  		ImageC<UIntT> lip_bin = GetBinaryLipImage(lip_pix,botfacehalf);
  		//Load up Ground Truth Image
  		ImageC<UIntT> gt_orig;
  		if(!Load(gtimg, gt_orig)) cerr<<"Could not load original Ground Truth Image"<<endl;
  		ImageC<UIntT> gt(gt_orig,botfacehalf);
  		RealT score = ComputeSegmentationQuality(lip_bin, gt);
  		cout<<score<<endl;
	}
	return 0;
}

template <typename T>
DListC<Index2dC> ClusterIdent(Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > &finalpop, const ImageC<T> &img, const IntT &segtype)
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

