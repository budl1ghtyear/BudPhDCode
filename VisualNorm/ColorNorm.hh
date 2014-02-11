#ifndef COLORNORM_HEADER
#define COLORNORM_HEADER 1

////////////////////////////////////////////////////////////////////
//! file="ColorNorm.hh"
//! author="Xuan Zou"
#include "Ravl/Point3d.hh"
#include 	"Ravl/Image/Image.hh"
#include 	"Ravl/Image/ByteRGBValue.hh"
#include 	"Ravl/Image/RealRGBValue.hh"
#include		"Ravl/Vector.hh"
#include 	"Ravl/Affine3d.hh"
#include		"Ravl/Array1dIter.hh"
#include 	"Ravl/Image/ImageRectangle.hh"


//! userlevel=Normal
//: Normalise the color of an image

using namespace RavlN;
using namespace RavlImageN;

enum ColorConstType
{  GREY_WORLD, // grey_world assumption, 
    CHANNEL_MAX,  //channel_ max assumption
    GREY_EDGE //grey_edge assumption
 };

//: Using Color Constancy method to estimate the color of illumination  and the reflectance of the object in the image
//: 
class ColorNormC
{	public: 
        	ColorNormC()
        	{};//: default constructor
 		
		VectorC IllumEstimate(const ImageC<ByteRGBValueC> &im, bool verbose );
                //: return estimated illumination based on color constancy method
		//!param: im - input color image

		VectorC IllumEstimate(const ImageC<ByteRGBValueC> &im, const ImageRectangleC &imRec , bool verbose );
		//: return the estimated illumination based on the img window
		//!param: im - input color image
		//!param: imRec - window of the interested region

		ImageC<ByteRGBValueC> Apply(const ImageC<ByteRGBValueC> &im, bool verbose );
		//: Apply color constancy method to image and obtain the corrected image
		//!param: im - input color image
		
		ImageC<ByteRGBValueC> Transform(const ImageC<ByteRGBValueC> &im, const VectorC illumSrc, const VectorC illumDes, bool verbose  );
		//: Apply color constancy method to image with the given illumination and obtain the transformed image
		//!param: im - input color image
		//!param: illumDes - desired illumination
		//!param: illumSrc - source illumination

		ImageC<ByteRGBValueC> Transform(const ImageC<ByteRGBValueC> &im, const VectorC illumDes, bool verbose  );
		//: Apply color constancy method to image with the given illumination and obtain the transformed image
		//!param: im - input color image
		//!param: illumDes - desired illumination

		bool SetColorConstType(ColorConstType type)
		{ m_type = type; return 1; }
                //: Set the type of color constancy method
                //!param: type - type of color constancy method

		ColorConstType GetColorConstType()
		{ return m_type; }
		//: Get the current type of color constancy method
		
	protected:  
		ColorConstType m_type;
		//: the current type of color constancy method

};


//: transform color based on affine transform in color space. The affine transform can be learnt from two exmaple images
//: 
class ColorNormAffineC
{	public: 
		ColorNormAffineC()
		{ m_affine = Affine3dC(); };
		//: Constructor

		bool SetTransform(Affine3dC& affine)
		{m_affine = affine; return 1; };
		//: Set affine transform
		//!param: affine - affine model

		Affine3dC GetTransform()
		{ return m_affine;}
		//: Get the current affine transform
               
		bool LearnTransformFromImgs(const ImageC<ByteRGBValueC>& imSrc, const ImageC<ByteRGBValueC>& imDes, bool verbose );
		//: learn affine transform from two example images and set the learnt transform to the model.
		//!param: imSrc - source color image
		//!param: imDes - target color image
		bool LearnTransformFromImgs(const ImageC<RealRGBValueC>& imSrc, const ImageC<RealRGBValueC>& imDes, bool verbose );
		

		Affine3dC LearnTransform(const ImageC<RealRGBValueC>& imSrc, const ImageC<RealRGBValueC>& imDes, bool verbose );
              	//: return the affine transform learnt from two explae images
		//!param: imSrc - source color image
		//!param: imDes - target color image

		ImageC<ByteRGBValueC> Apply(const ImageC<ByteRGBValueC>& im, bool verbose );
		//: return the affine transformed color image
		//!param: im - input color image
		ImageC<RealRGBValueC> Apply(const ImageC<RealRGBValueC>& im, bool verbose );

	protected:
		Affine3dC m_affine;
               //: affine transform in color space
};


//: learn a color transform so that the transformed source image wil have same mean and std with the destination image
class ColorNormHistStatC
{	public: 
		ColorNormHistStatC()
		{}//: constructor	

		bool GetHistStatistics(const ImageC<ByteRGBValueC> &im, Vector3dC &meanVec,  Vector3dC &stdVec);
		//: get the statistics for the image 
		//!param: im - input image
		//!param: meanVec - mean vector (R,G,B) of the img;
		//!param: stdVec - std vector of the img;

		bool GetHistStatistics(const ImageC<ByteRGBValueC> &im, const ImageRectangleC &imRec, Vector3dC &meanVec,  Vector3dC &stdVec);
		//: get the statistics of a specfic region of an image
		//!param: im - input image
		//!param: imRec - interested region
		//!param: meanVec - mean vector (R,G,B) of the img;
		//!param: stdVec - std vector of the img;

		ImageC<ByteRGBValueC> Apply(const ImageC<ByteRGBValueC> &imRGB, const Vector3dC &meanVecDes, const Vector3dC & stdVecDes, const Vector3dC &meanVecSrc, const Vector3dC& stdVecSrc);
		//: learn the transform based on the mean and std vectors of the source and destination image, and apply it to the input image 
		//!param: imRGB - input image, 
		//!param: meanVecDes - mean vec of destination image, 
		//!param: stdVecDes - std vec of the destination image,
		//!param: meanVecSrc - mean vec of the source image,
		//!param: stdVecSrc - std vec of the source image
	
};

//: learn a color transform so that the transformed source image wil have same histogram with the destination image
class ColorNormHistShapeC
{	public:
		ColorNormHistShapeC()
		{}; //: constuctor

		bool GetColorMap(const ImageC<ByteT> &imSrc, const ImageC<ByteT> &imDes, Array1dC<ByteT> &m_ColorMap);		
		//: learn the color mapping from the source and destination images
		//!param: imSrc - source image
		//!param: imDes - destination image
		//!param: m_ColorMap - color mapping learnt 

		bool GetColorMap(const ImageC<ByteRGBValueC> &imSrc, const ImageC<ByteRGBValueC> &imDes);
		//: learn the color mapping from the source and destination images
		//!param: imSrc - source image
		//!param: imDes - destination image
		
		bool GetColorMap(const ImageC<ByteRGBValueC>&imSrc,const ImageRectangleC &imRecSrc, const ImageC<ByteRGBValueC> &imDes, const ImageRectangleC &imRecDes);
		//: learn the color mapping from the corresponding regions in source and destination images
		//!param: imSrc - source image
		//!param: imRecSrc - interested region in source image
		//!param: imDes - destination image
		//!param: imRecDes - interested region in destination image

		ImageC<ByteRGBValueC> Apply(const ImageC<ByteRGBValueC> &imRGB);
		//: Apply the learnt color mapping to the input image
		//!param: imRGB - input image
		
	protected:
		Array1dC<ByteT> m_colorMapR; //:color mapping array for R channel
		Array1dC<ByteT> m_colorMapG; //:color mapping array for G channel
		Array1dC<ByteT> m_colorMapB; //:color mapping array for B channel
};


#endif
