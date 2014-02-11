#ifndef LIPSEGC_HH
#define LIPSEGC_HH
//////////////////////////////////////////////////////////////////////////////
//! file =      "Ravl/LipCode/Segmentation/SegmentationLibrary/LipSegmentation.hh"   
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
//BUDSW
#include "MouthRegion.hh"
#include "ColourConvert.hh"
#include "CostFunctions.hh"
#include "LipClusteringC.hh"
#include "ColourSpaceImage.hh"
using namespace RavlN;
using namespace RavlImageN;


//: Class to encapsulate the process of lip segmentation as detailed in B. Goswami lip segmentation publications
class LipSegmentationC
{
	public:
	LipSegmentationC(){}
	~LipSegmentationC(){}
	LipSegmentationC(const FilenameC &file);
	//:Constructor for inputs as file names
	LipSegmentationC(const ImageC<RealRGBValueC> &img);
	//:Constructor for inputs as images (must always use ImageC<RealRGBValueC>)
	
	//Apply various stages
	ColourSpaceImage ColourSpaceConversion(const PixelType &ptype);
	//: First stage is to get the colour conversion performed depending on specified pixel type, returns a ColourSpaceImage as its output
	Tuple4C<Affine2dC, ImageRectangleC, ImageRectangleC,ImageC<VectorC> > AffineMouthRegionDetection(const ColourSpaceImage &im);
	//: This method performs Affine normalisation of the face image and eye detection. This is used together to perform mouth-region detection. The first ImageRectangle is the mouth-region in affine space and the second is in normal image space
	Tuple2C<MeanCovarianceC,MeanCovarianceC> LipClusteringStep(const ClusteringType &type, const ImageC<VectorC> &im,const RealT &hval, const UIntT &nclust);
	//: This method performs clustering depending on a specified clustering method. See LipClusteringC for more....
	ImageC<UIntT>  ProbabilisticPixelLabellingStep(const ImageC<VectorC> &im, const MeanCovarianceC &sk, const MeanCovarianceC &lp, const LabellingType &type);
	//: ProbabilisticPixelLabellingStep() will label an image using some log-likelihood function provided the MeanCovariance estimates of the 2 classes are provided
	ImageC<UIntT> RegionIdentificationStep(const ImageC<UIntT> &im, const RegionIdentificationType &type);
	//: RegionIdentificationStep() this will perform lip region identification following pixel labelling using a specified region identification algorithm
	Tuple3C<ImageC<UIntT>, MeanCovarianceC, MeanCovarianceC> Apply(const PixelType &ptype,const ClusteringType &ctype,const LabellingType &ltype,const RegionIdentificationType &rtype, const RealT &hval, const UIntT &nclust);
	
	//:Accessor Method
	ImageC<RealRGBValueC> GetOriginalImage(void) const {return inp_img;}
	Affine2dC GetAffineTransform(void) const {return aff_transform;}
	Tuple2C<MeanCovarianceC,MeanCovarianceC> GetColorTrends(void) const {return color_trends;}
	protected:
	ImageC<RealRGBValueC> inp_img;	
	FilenameC inp_img_name;
	Affine2dC aff_transform;
	Tuple2C<MeanCovarianceC,MeanCovarianceC> color_trends;
};
#endif
