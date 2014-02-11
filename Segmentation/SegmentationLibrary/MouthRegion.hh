#ifndef MouthRegion_HH
#define MouthRegion_HH
//////////////////////////////////////////////////////////////////////////////
//! file =      "Ravl/LipCode/Segmentation/SegmentationLibrary/MouthRegion.hh"   
//! author =    "Bud Goswami"
//! lib =       SegLib
//! date =      "Ongoing"
//! userlevel = Normal
//! docentry =  "Ravl.API.LipCode.Segmentation.SegmentationLibrary"
#include "Ravl/Affine2d.hh"
#include "Ravl/Array1d.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImageConv.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/DrawFrame.hh"
#include "Ravl/Image/GaussConvolve2d.hh"
#include "Ravl/Image/OpenCVConvert.hh"
#include "Ravl/Image/ImageRectangle.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Image/WarpAffine.hh"
#include "Ravl/IO.hh"
#include "Ravl/PatternRec/Sample.hh"
#include "Ravl/Option.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/String.hh"
#include "Ravl/Tuple2.hh"
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

using namespace RavlImageN;
using namespace RavlConstN;
using namespace OmniN;
//:Methods related to Mouth Region Extraction


void detect_and_draw( IplImage* img, CvHaarClassifierCascade* cascade  );
//:Function prototype for detecting and drawing an object from an image using OpenCV
ImageRectangleC GetFaceCoords(const FilenameC &img_name, const FilenameC &c_name);
//:Wrapper Function to use OpenCv to get Facerectangle
ImageRectangleC GetFaceCoords(const ImageC<RealRGBValueC> &img, const FilenameC &c_name);
//:Same Function but using an Image as an input. PLEASE NOTE that the RAVL OpenCV Image converter at the time of documentation was buggy. Don't use this method!
//OMNI Functions
Tuple2C<Point2dC,Point2dC> GetEyeLoc(const FilenameC &ImageFile);
//:Omni SW to obtain eye centre positions
Tuple2C<Point2dC,Point2dC> GetEyeLoc(const ImageC<RealRGBValueC> &img);
//:Same Function but using Image as an input
//NORMAL RAVL STUFF
ImageRectangleC GetMouthRegion(const FilenameC &img_name, const FilenameC &c_name, const IndexC &window);
//:GetMouthRegion() using image filename, the cascade name and a padding window
ImageRectangleC GetMouthRegion(const ImageC<RealRGBValueC> &img, const FilenameC &c_name, const IndexC &window);
//:GetMouthRegion() using image, the cascade name and a padding window
template <typename T>
Tuple2C<ImageC<T>,Affine2dC> AffNormWRTEyes(const ImageC<T> & image, const Point2dC &nle,const Point2dC &nre)
//:Affine normalise an image wrt to detected eye positions
{
    ImageC<T> normalisedImage;
/*
//!Gabor
    RealT rows = 160;
    RealT cols = 128;
//for leye=28c,55r
    RealT rowFrac = 0.34375;
    RealT colFrac = 0.21875;
*/
  RealT rows = 142;
  RealT cols =  120;
  RealT rowFrac = 0.35;
  RealT colFrac = 0.25;

//for FRGC 3d Data
//        RealT rowFrac = 0.3;
 //   RealT colFrac = 0.1;
    UIntT gaussOrder = 7;
        
    //: Workout a simple geometric normalisation
    //==========================================
    Point2dC le(rows * rowFrac, (1.0-colFrac) * cols);
    Point2dC re(rows * rowFrac, cols * colFrac);
    ImageRectangleC outRect(Round(rows), Round(cols));
    Vector2dC d1 = le - re;
    Vector2dC d2 = nle - nre;
    RealT rot = d2.Angle() - d1.Angle();     
    RealT scale = d2.Modulus() / d1.Modulus();  
    Matrix2dC rotm = Matrix2dC(Cos(rot) * scale,-Sin(rot) * scale, Sin(rot) * scale,Cos(rot) * scale);      
    Point2dC cent= ((nle + nre)/2);
    Point2dC dcent = rotm * (le + re)/2;
    Vector2dC off = (cent - dcent);
    Affine2dC tr(rotm,off);
    
    //: Smooth and geometric normalise
    //=================================
    //GaussConvolve2dC<T> smooth(gaussOrder);
    WarpAffineC<T> warp(outRect,tr);
    ImageRectangleC irec = warp.InputRectangle();
    irec = irec.Expand(gaussOrder/2 + 2);
    irec.ClipBy(image.Frame());
    //normalisedImage = smooth.Apply(ImageC<T>(image,irec));
	normalisedImage = ImageC<T>(image,irec).Copy();
//  warp.SetOutputRectangle(IndexRange2dC(144,120));
    normalisedImage = warp.Apply(normalisedImage);
    
    //: Lets display the image
    //========================
	Tuple2C<ImageC<T>,Affine2dC> out(normalisedImage, tr);
   return out;
}
//:Affine Normalise the image wrt to standard eye positions

template <typename T>
SampleC<VectorC> GetSample(const ImageC<T> &img)
//: Convert an image to a SampleC<VectorC>
{
	SArray1dC<VectorC> sar(img.Frame().Area());
	UIntT dim = img[img.Frame().TopLeft()].Size();
	UIntT ind = 0;
	for(Array2dIterC<T> it(img); it; it++)
	{
		VectorC v(dim);
		for(UIntT i = 0; i < dim; i++)
		{
			v[i] = (*it)[i];
		}		
		sar[ind] = v.Copy();
		ind++;
	}
	SampleC<VectorC> out(sar);
	return out;
}


#endif

