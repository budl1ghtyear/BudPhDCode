#ifndef ProcrustesAlign_HH
#define ProcrustesAlign_HH

//INCLUDE STATEMENTS
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

class ProcrustesAlignC
{
	public:
	ProcrustesAlignC(){}
	~ProcrustesAlignC(){}
	//Constructor with the mean shape to align the shapes to:
	ProcrustesAlignC(const Array1dC<Point2dC> &ref);
	ProcrustesAlignC(const SArray1dC<Point2dC> &ref);
	ProcrustesAlignC(const VectorC &ref);
	//Apply method:
	DListC<SArray1dC<Point2dC> > Apply(const DListC<SArray1dC<Point2dC> > &pts);
	SArray1dC<Point2dC> AlignShape(const SArray1dC<Point2dC> &pts);
	//Accessor methods:
	SArray1dC<Point2dC> GetReferenceShape(void) const {return ref_shape;}
		
	protected:
	SArray1dC<Point2dC> ref_shape;
	SArray1dC<Point2dC> AffineProject(const SArray1dC<Point2dC> &pts, const Affine2dC &aff);
	SArray1dC<Point2dC> ComputeMeanShape(const DListC<SArray1dC<Point2dC> > &pts);
};


#endif
