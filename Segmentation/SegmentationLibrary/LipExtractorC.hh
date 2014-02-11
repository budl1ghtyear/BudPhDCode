#ifndef LIPEXTRACTORC_HH
#define LIPEXTRACTORC_HH
//////////////////////////////////////////////////////////////////////////////
//! file =      "Ravl/LipCode/Segmentation/SegmentationLibrary/LipExtractor.hh"   
//! author =    "Bud Goswami"
//! lib =       SegLib
//! date =      "10.08.09"
//! userlevel = Normal
//! docentry =  "Ravl.API.LipCode.Segmentation.SegmentationLibrary"


//INCLUDE FILES
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/SArray1d.hh"
#include "Ravl/SArray1dIter3.hh"
#include "Ravl/Option.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/String.hh"
#include "LipSegmentation.hh"
#include "ColourSpaceImage.hh"
#include "XMLToGTConv.hh"
#include "BANCAGT.hh"
#include "Ravl/Boundary.hh"
#include "Ravl/Image/MorphClose.hh"
#include "Ravl/Image/Erode.hh"
#include "Ravl/Ellipse2d.hh"
#include "Ravl/Image/DrawEllipse.hh"
#include "Ravl/Array2dIter3.hh"
#include "Ravl/Math.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/LineABC2d.hh"
#include "Ravl/Complex.hh"
#include "Ravl/Image/DrawLine.hh"
#include "Ravl/Pair.hh"
#include "Ravl/Image/ImageRectangle.hh"
#include "Ravl/Boundary.hh"
using namespace RavlN;
using namespace RavlImageN;
using namespace RavlConstN;
//:Purpose - to get the lip shape boundary from a binary lip boundary image
class LipExtractorC
{
	//PUBLIC
	public:
	LipExtractorC() {}
	~LipExtractorC() {}
	LipExtractorC(const ImageC<UIntT> &img, const ImageRectangleC &mR, const UIntT &num_pts);
	//: Constructor with the image, the ROI and the number of points extracted <p>
	//: Accessor Methods - 
	ImageC<UIntT> GetOriginalLipBoundaryImage(void) const {return lip_img;}
	
	UIntT GetNumPoints(void) const {return num_points;}
	Ellipse2dC GetLipEllipse(void) const {return lip_ellipse;}
	SArray1dC<Point2dC> GetAllBoundaryPoints(void) const {return lip_bound_pts;}
	SArray1dC<Point2dC> GetLipBoundaryPoints(void) {return selected_Points;}
	Array1dC<Point2dC> GetLipBoundaryPointsArray(void) { return Array1dC<Point2dC>(selected_Points);}
	BoundaryC GetOuterLipBoundaryPixelList(void) const {return outerLipBoundary;}
	ImageC<UIntT> PerformMorphology(void);
	//PROTECTED VARIABLES AND DATA
	protected:
	ImageC<UIntT> lip_img;
	UIntT num_points; //:THIS NEEDS TO BE AN EVEN NUMBER BECAUSE OF OUR IMPLEMENTATION
	UIntT NUM_POINTS_EACH_SIDE; //:Automatically calculate inside the constructor
	SArray1dC<Point2dC> selected_Points; //:For the final selected lip boundary points
	SArray1dC<Point2dC> lip_bound_pts; //:Contains all the lip boundary points
	Ellipse2dC lip_ellipse;//: Best fitting ellipse for our lip points
	SArray1dC<Vector2dC> allVectors;
	PairC<IndexC> lip_corner_points;
	ImageRectangleC mouthRegion;
	BoundaryC outerLipBoundary;
	
	SArray1dC<Point2dC> Apply(void);
	
	PairC<IndexC> GetLipCornerPoints(void);
	SArray1dC<Vector2dC> GetLipNormalVectors(void);
	SArray1dC<Point2dC> SelectLipBoundaryPoints(void);
};


#endif
