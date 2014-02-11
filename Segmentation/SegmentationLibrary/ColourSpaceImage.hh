#ifndef COLOURSPACEIMAGE_HH
#define COLOURSPACEIMAGE_HH
//////////////////////////////////////////////////////////////////////////////
//! file =      "Ravl/LipCode/Segmentation/SegmentationLibrary/ColourSpaceImage.hh"   
//! author =    "Bud Goswami"
//! lib =       SegLib
//! date =      "10.08.09"
//! userlevel = Normal
//! docentry =  "Ravl.API.LipCode.Segmentation.SegmentationLibrary"
#include "Ravl/Array1d.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"  
#include "Ravl/Image/ConnectedComponents.hh"
#include "Ravl/Image/OpenCVConvert.hh"
#include "Ravl/Image/ImageRectangle.hh"
#include "Ravl/Image/Segmentation.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/SumsNd2.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImageRectangle.hh"
#include "Ravl/Array2dIter2.hh"
#include "Ravl/Vector.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/Image/RealHSVValue.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Image/RealYUVValue.hh"
#include "Ravl/Tuple2.hh"
#include "Ravl/Tuple3.hh"
#include "Ravl/Tuple4.hh"

typedef enum {RGB=0, NormRG=1, HSV=2, YUV=3, CIE=4, PHue=5} PixelType;	
//: Enumerated type for pixel type
typedef enum {JNormal=0, JSkin=1} LabellingType;
//: Enumerated type for Pixel Labelling type
typedef enum {CC=0,CCentral=1} RegionIdentificationType;
//; Enumerated type for Region Identification Type

//: Encapsulating data structure required by the lip segmentation algorithm
class ColourSpaceImage
{
	public:
	ColourSpaceImage(){}
	~ColourSpaceImage(){}
	template <class T>
	ColourSpaceImage(const ImageC<T> &im, const PixelType &c_sp)
	//: Constructor with image and the pixel type
	{
		dim = im[im.Frame().TopLeft()].Size();
		ImageC<VectorC> vec_im(im.Frame());
		for(Array2dIter2C<VectorC,T> it(vec_im,im); it; it++)
		{
			VectorC vec(dim);
			for(UIntT i = 0; i < dim; i++)
			{
				vec[i] = it.Data2()[i];
			}
			it.Data1() = vec.Copy();
		}
		img = vec_im.Copy();
		c_space = c_sp;
	}
	//Accessor Methods
	ImageC<VectorC> GetImage() const {return img;}
	PixelType GetImageType() const {return c_space;}
	UIntT GetNumChannels() const {return dim;}	
	protected:
	ImageC<VectorC> img; //:to store the image that we have
	PixelType c_space; //:to store the colour space
	UIntT dim; //:number of colour channels
};

#endif
