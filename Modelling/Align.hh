#ifndef ALIGN_HH
#define ALIGN_HH
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/SArray1d.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/Vector.hh"
#include "Ravl/Affine2d.hh"
#include "Ravl/Image/Image.hh"
using namespace RavlN;
using namespace RavlImageN;


//: Class to perform single pass similarity or affine alignment
class AlignC
{
	public:
	AlignC(){}
	~AlignC(){}
	//Apply methods:
	DListC<SArray1dC<Point2dC> > PerformSimilarityAlignment(const DListC<SArray1dC<Point2dC> > &pts, const SArray1dC<Point2dC> &ref);
	DListC<SArray1dC<Point2dC> > PerformAffineAlignment(const DListC<SArray1dC<Point2dC> > &pts, const SArray1dC<Point2dC> &ref);
	DListC<SArray1dC<Point2dC> > ComputeDeviation(const DListC<SArray1dC<Point2dC> > &pts,  const SArray1dC<Point2dC> &ref);
	protected:
	SArray1dC<Point2dC> AffineProject(const SArray1dC<Point2dC> &pts, const Affine2dC &aff);
};

#endif
