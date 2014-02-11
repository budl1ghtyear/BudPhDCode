#ifndef BANCAGT_HH
#define BANCAGT_HH
/*
 * Class - BANCAGT.hh
 * Author - Bud Goswami
 * Date - 10.08.09
 */
//////////////////////////////////////////////////////////////////////////////
//! file =      "Ravl/LipCode/Segmentation/SegmentationLibrary/BANCAGT.hh"   
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
#include "Ravl/Stream.hh"
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------
// **********  BANCAGT  **********************************************
// ---------------------------------------------------------------------------
//:Class to generate ground truth data from a BANCA annotation (.did) file <p>
//: Note that the BANCA .did files are in a strange order and this class compensates for that
class BANCAGT
{
	public:
	BANCAGT(){}
	~BANCAGT(){}
	
	BANCAGT(const FilenameC &fn, const ImageRectangleC &imrec);
	//: Constructor with filename and an imagerectangle (to specify region of interest)
	//Accessor Methods
	StringC GetImageNumber(void) const {return f_name.BaseNameComponent();}
	ImageC<UIntT> GetGTImage(void) const {return lip_bin;}
	ImageC<UIntT> Apply(void); 
	DListC<Point2dC> GetOuterPoints(void) const {return this->outer;}
	DListC<Point2dC> GetInnerPoints(void) const {return this->inner;}
	protected:
	FilenameC f_name; //: Store the filename of the original image
	DListC<Point2dC> outer,inner; //: Store the outer and inner lip co-ordinates as obtained from the .did file
	ImageC<UIntT> lip_bin; //:to store the binary lip image
	DListC<Point2dC> RearrangePoints(DListC<Point2dC> &pts);
	//: The .did files are in an order where all the x co-ordinates are listed first before the y co-ordinates.<p>
	//: This method reorders the data in the .did file to present meaningful and correct co-ordinates
};
#endif
