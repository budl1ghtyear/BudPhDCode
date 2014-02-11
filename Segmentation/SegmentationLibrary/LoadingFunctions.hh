//ImageLoadingFunctions.hh
//Author - Bud Goswami
//Purpose - to create a header file that contains the function definitions of colour conversion using RAVL
#ifndef LoadingFunctions_HH
#define LoadingFunctions_HH

//DATA STRUCTURES
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter2.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Vector.hh"
#include "Ravl/Vector2d.hh"
//IMAGE STUFF
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/OpenCVConvert.hh"
#include "Ravl/Image/ImageRectangle.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/IO.hh"
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
#include "CascadedMCD.hh"
#include "KSmirnovMCD.hh"
using namespace RavlImageN;
using namespace RavlConstN;

ImageC<RealRGBValueC > LoadImage(const ImageC<RealRGBValueC> &imfile, const UIntT &c_opt, const ImageRectangleC &botfacehalf)
{
		ImageC<RealRGBValueC> m_img(GetRGBVals(imfile),botfacehalf);
		return m_img;
 
}

ImageC<RealRGBValueC > LoadImage(const FilenameC &imfile, const UIntT &c_opt, const ImageRectangleC &botfacehalf)
{
		ImageC<RealRGBValueC> m_img(GetRGBVals(imfile),botfacehalf);
		return m_img;
}

templat


/*
//Using the provided file name, load up the correct image that we want
template <typename T>
ImageC<T> LoadImage(const ImageC<RealRGBValueC> &imfile, const UIntT &c_opt, const ImageRectangleC &botfacehalf)
{
	switch(c_opt)
  	{
  		case 0:
  		{
  			ImageC<RealRGBValueC> m_img(GetRGBVals(imfile),botfacehalf);
  			return m_img;
  			break;
  		}
 		case 1:
  		{
  			ImageC<VectorC> m_img(GetNormRGVals(imfile),botfacehalf);
			return m_img;
  			break;
  		}
  		case 2:
  		{
  			ImageC<RealHSVValueC> m_img(GetHSVVals(imfile),botfacehalf);
			return m_img;
			break;
  		}
  		case 3:
  		{
  			ImageC<RealYUVValueC> m_img(GetYUVVals(imfile),botfacehalf);
			return m_img;
   		break;
  		}
 		case 4:
  		{
  			ImageC<VectorC> m_img(GetCIELabVals(imfile),botfacehalf);
			return m_img;				
			break;
		}
		case 5:
		{
			ImageC<VectorC> m_img(GetPseudoHueVals(imfile),botfacehalf);
			return m_img;
			break;
		}
		default:
		{
			ImageC<RealRGBValueC> m_img(GetRGBVals(imfile),botfacehalf);
  			return m_img;
  			break;		
		}
	}	
}


//Using the provided file name, load up the correct image that we want
template <typename T>
ImageC<T> LoadImage(const FilenameC &imfile, const UIntT &c_opt, const ImageRectangleC &botfacehalf)
{
	switch(c_opt)
  	{
  		case 0:
  		{
  			ImageC<RealRGBValueC> m_img(GetRGBVals(imfile),botfacehalf);
  			return m_img;
  			break;
  		}
 		case 1:
  		{
  			ImageC<VectorC> m_img(GetNormRGVals(imfile),botfacehalf);
			return m_img;
  			break;
  		}
  		case 2:
  		{
  			ImageC<RealHSVValueC> m_img(GetHSVVals(imfile),botfacehalf);
			return m_img;
			break;
  		}
  		case 3:
  		{
  			ImageC<RealYUVValueC> m_img(GetYUVVals(imfile),botfacehalf);
			return m_img;
   		break;
  		}
 		case 4:
  		{
  			ImageC<VectorC> m_img(GetCIELabVals(imfile),botfacehalf);
			return m_img;				
			break;
		}
		case 5:
		{
			ImageC<VectorC> m_img(GetPseudoHueVals(imfile),botfacehalf);
			return m_img;
			break;
		}
		default:
		{
			ImageC<RealRGBValueC> m_img(GetRGBVals(imfile),botfacehalf);
  			return m_img;
  			break;		
		}
	}	
}
*/
//Method if we want our original slow output
DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > MCDClustering(const FilenameC &im, const UIntT &c_type, const UIntT &im_type, const ImageRectangleC &irec);

//Method for fast output
DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > MCDClustering(const ImageC<RealRGBValueC> &im, const UIntT &c_type, const UIntT &im_type, const ImageRectangleC &irec);
#endif

