#ifndef SEGC_HH
#define SEGC_HH

/*
 * CLASS - SegmentationC & ColourSpaceImage
 * Author - Bud Goswami
 * Date - 10.08.09 
 */
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"  
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImageRectangle.hh"
#include "Ravl/Array2dIter2.hh"
#include "Ravl/Vector.hh"
#include "Ravl/Image/RealHSVValue.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Image/RealYUVValue.hh"
#include "Ravl/Tuple4.hh"
//BUDSW
#include "MouthRegion.hh"
#include "ColourConvert.hh"
#include "GlobalVariables.hh"

using namespace RavlN;
using namespace RavlImageN;
class ColourSpaceImage
{
	public:
	ColourSpaceImage(){}
	~ColourSpaceImage(){}
	template <class T>
	ColourSpaceImage(const ImageC<T> &im, ImType &c_sp)
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
	enum ImType{RGB,NRG,HSV,YUV,CIELAB,PHUE};
	//Accessor Methods
	ImageC<VectorC> GetImage() const {return img;}
	ImType GetImageType() const {return c_space;}
	UIntT GetNumChannels() const {return dim;}	
	protected:
	ImageC<VectorC> img; //to store the image that we have
	ImType c_space; //to store the colour space
	UIntT dim; //number of colour channels
};

class SegmentationC
{
	public:
	//ENUMERATED TYPES FOR VARIOUS OPTIONS
	enum PixelType {RGB, NormRG, HSV, YUV, CIE, PHue};
	enum ClusteringType {CMCD, KSMCD, KMeans, FCM};
	enum LabellingType {JNormal, JSkin};
	enum RegionIdenticationType {CC}; //Connected Components
	
	SegmentationC(){}
	~SegmentationC(){}
	SegmentationC (const FilenameC &file, PixelType &ptype, ClusteringType &ctype, LabellingType &ltype, RegionIdenticationType &rtype, const RealT &hvalue);
	SegmentationC (const FilenameC &file, PixelType &ptype, ClusteringType &ctype, LabellingType &ltype, RegionIdenticationType &rtype);
	Tuple3C<ImageC<UIntT>, MeanCovarianceC, MeanCovarianceC> Apply(void);
	//Accessor Methods	
	ColourSpaceImage GetColourSpaceImage() const {return my_im;}
	MeanCovarianceC GetLipColourEstimate() const {return lip;}
	MeanCovarianceC GetSkinColourEstimate() const {return skin;}
	ImageC<UIntT> GetLabelledImage() const {return labelled_im;}
	ImageC<UIntT> GetLipBinaryImage() const {return lipbinary_im;}
	protected:
	FilenameC in_file;
	ImageC<RealRGBValueC> input_img;
	ColourSpaceImage my_im;
	MeanCovarianceC lip;
	MeanCovarianceC skin;
	ImageC<UIntT> labelled_im;
	ImageC<UIntT> lipbinary_im;
	ImageRectangleC origmouth_reg;
	ImageRectangleC affmouth_reg;
	Affine2dC aff;
	//PROTECTED METHODS - Implement the modules of lip segmentation
	Tuple3C<Affine2dC, ImageRectangleC, ImageRectangleC> AffineMouthRegionDetection(const ImageC<RealRGBValueC> &im, const PixelType &ptype);
	//:Get Affine Normalised Image For Use In Clustering
	Tuple2C<MeanCovarianceC, MeanCovarianceC> GetChromaticClusters(ColourSpaceImage &im, const ClusteringType &ctype, const RealT &hval);
	//:Get Cluster Trends Info For Probabilistic Pixel Labelling
	ImageC<UIntT> ProbabilisticPixelLabellingStep(const ImageC<T> &im, const MeanCovarianceC &sk, const MeanCovarianceC &lp, const UIntT &type);
	//:Perform Probabilistic Pixel Labelling
	ImageC<UIntT> RegionIdenticationStep(const ImageC<UIntT> &im, const ImageRectangleC &imrec, const RegionIdenticationType &rtype);
	//:Perform Region Identification
	//Utility Methods
};
#endif
