#ifndef XMLCONVC_HH
#define XMLCONVC_HH
//////////////////////////////////////////////////////////////////////////////
//! file =      "Ravl/LipCode/Segmentation/SegmentationLibrary/XMLToGTConv.hh"   
//! author =    "Bud Goswami"
//! lib =       SegLib
//! date =      "10.08.09"
//! userlevel = Normal
//! docentry =  "Ravl.API.LipCode.Segmentation.SegmentationLibrary"
#include "Ravl/Option.hh"
#include "Ravl/Image/ImagePointFeatureSet.hh"
#include "Ravl/HashIter.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/StrStream.hh"
#include "Ravl/XMLStream.hh"
#include "Ravl/EntryPnt.hh"
#include "Ravl/SArray1d.hh"
#include "Ravl/SArray1dIter.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/Index2d.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/DP/FileFormatIO.hh"
#include "Ravl/Polygon2d.hh"

using namespace RavlN;
using namespace RavlImageN;
//: Class to generate Ground Truth image given some XML GT files
class XMLToGTConv
{
	public:
	XMLToGTConv(){}
	~XMLToGTConv(){}
	
	XMLToGTConv(const FilenameC &fn, const ImageRectangleC &imrec);
	//: Constructor with XML ground truth filename and the ROI in the original image
	//Accessor Methods
	StringC GetImageNumber(void) const {return f_name.BaseNameComponent();}
	ImagePointFeatureSetC GetImageFeaturePointSet(void) const {return im_data;}
	ImageC<UIntT> GetGTImage(void) const {return lip_bin;}
	ImageC<UIntT> Apply(void); 
	//:Apply method
	protected:
	FilenameC f_name;
	ImagePointFeatureSetC im_data; //:to store the xml contents
	ImageC<UIntT> lip_bin; //:to store the binary lip image
};
#endif
