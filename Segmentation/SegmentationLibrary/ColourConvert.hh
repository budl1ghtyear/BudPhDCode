//ColourConvert.hh
//Author - Bud Goswami
//Purpose - to create a header file that contains the function definitions of colour conversion using RAVL
#ifndef ColourConvert_HH
#define ColourConvert_HH
//////////////////////////////////////////////////////////////////////////////
//! file =      "Ravl/LipCode/Segmentation/SegmentationLibrary/ColourConvert.hh"   
//! author =    "Bud Goswami"
//! lib =       SegLib
//! date =      "10.09.09"
//! userlevel = Normal
//! docentry =  "Ravl.API.LipCode.Segmentation.SegmentationLibrary"
//DATA STRUCTURES
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter2.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/Vector.hh"
#include "Ravl/Vector2d.hh"
//IMAGE STUFF
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/ImageConv.hh"
#include "Ravl/Image/RealHSVValue.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Image/RealYUVValue.hh"
#include "Ravl/IO.hh"
#include "Ravl/Option.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/StdMath.hh"
#include "Ravl/Matrix.hh"
#include "Ravl/Math.hh"

using namespace RavlImageN;
using namespace RavlConstN;
//: Methods to convert between pixel types from RGB. Method names are self-explanatory
ImageC<RealRGBValueC> GetRGBVals(const FilenameC &imgname);
ImageC<RealRGBValueC> GetRGBVals(const ImageC<RealRGBValueC> &img);
ImageC<VectorC> GetNormRGVals(const FilenameC &imgname);
ImageC<VectorC> GetNormRGVals(const ImageC<RealRGBValueC> &img);
ImageC<RealHSVValueC> GetHSVVals(const FilenameC &imgname);
ImageC<RealHSVValueC> GetHSVVals(const ImageC<RealRGBValueC> &img);
ImageC<RealYUVValueC> GetYUVVals(const FilenameC &imgname);
ImageC<RealYUVValueC> GetYUVVals(const ImageC<RealRGBValueC> &img);
ImageC<VectorC> GetCIELabVals(const FilenameC &imgname);
ImageC<VectorC> GetCIELabVals(const ImageC<RealRGBValueC> &img);
ImageC<VectorC> GetPseudoHueVals(const FilenameC &imgname);
ImageC<VectorC> GetPseudoHueVals(const ImageC<RealRGBValueC> &img);
template <typename T>
ImageC<VectorC> ConvertToVecImage(const ImageC<T> &im)
//: Pixel structure converter, ImageC<T> -> ImageC<VectorC>
{
	UIntT dim = im[im.Frame().TopLeft()].Size();
	ImageC<VectorC> vec_im(im.Frame());
	for(Array2dIter2C<T, VectorC> it(im, vec_im); it; it++)
	{
		VectorC vec(dim);
		for(UIntT i = 0; i < dim; i++)
		{
			vec[i] = it.Data1()[i];		
		}
		it.Data2() = vec.Copy();
	}
	return vec_im;
}
#endif
