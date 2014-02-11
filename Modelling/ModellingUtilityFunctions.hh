#ifndef MODELLINGUTILITYFUNCTIONS_HH
#define MODELLINGUTILITYFUNCTIONS_HH

//INCLUDE FILES
//Required Libs
#include "Ravl/Affine2d.hh"
#include "Ravl/Array2d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/SArray2d.hh"
#include "Ravl/SArray1dIter.hh"
#include "Ravl/DArray1d.hh"
#include "Ravl/DP/SequenceIO.hh"
#include "Ravl/DP/FileFormatIO.hh"
#include "Ravl/IO.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/DrawCircle.hh"
#include "Ravl/Image/DrawLine.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/LinePP2d.hh"
#include "Ravl/Math.hh"
#include "Ravl/Matrix.hh"
#include "Ravl/Option.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/SArray1d.hh"
#include "Ravl/SArray2d.hh"
#include "Ravl/SArray1dIter.hh"
#include "Ravl/Stream.hh"
#include "Ravl/Tuple2.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Vector2d.hh"
#include "Ravl/Vector.hh"
#include "Ravl/Image/WarpAffine.hh"

using namespace RavlN;
using namespace RavlImageN;

class ModellingUtilityFunctions
{
	public:
	ModellingUtilityFunctions(){;}
	~ModellingUtilityFunctions(){;}
	
	//True Utility Functions to load data
	DListC<SArray1dC<Point2dC> > LoadTRESData(const FilenameC& f); //: Load data that has been produced with the TRES format
	ImageC<RealRGBValueC> LoadImage(const FilenameC& f); //: Load an image into an ImageC<RealRGBValueC>
	//Data-structure manipulation methods
	SArray1dC<Point2dC> ReorderPoints(const SArray1dC<Point2dC> &lp_pts); //:Reorder Spline Points for EJ Landmarked Data
	Array1dC<Point2dC> SArrayToArray(const SArray1dC<Point2dC> &pts); //: Convert between SArray1dC and Array1dC
	SArray1dC<Point2dC> ArrayToSArray(const Array1dC<Point2dC> &pts);
	VectorC ArrayToVector(const Array1dC<Point2dC> &arr); //: Convert an array to a vector
	//Normalisation Functions
	Vector2dC ScaleNormalise(const RealT &mn, const RealT &md);//: Perform scale normalisation
	RealT RotationNormalise(const Array1dC<Point2dC> &arr); //: Perform Rotation + Translation Normalisation
	Point2dC TranslationNormalise(const Array1dC<Point2dC> &arr); //: Perform Translation Normalisation
	Array1dC<Point2dC> AffineNormalise(const Array1dC<Point2dC> &pts, const Point2dC &trans, const Vector2dC &scale,const RealT &rotation);//:Perform some affine normalisation using supplied affine parameters
	//Data projection and statistics computing functions
	Array1dC<Point2dC> PerformProjection(const Affine2dC &aff, const Array1dC<Point2dC> &pts); //: Perform an affine projection using some point and an affine transform
	Point2dC ComputeMeanPoint(const Array1dC<Point2dC> &pts); //:Compute the mean point of a set of points
};
#endif
