#ifndef BaseStateVector_HH
#define BaseStateVector_HH


#include "Ravl/SArray1d.hh"
#include "Ravl/SArray1dIter.hh"
#include "Ravl/Vector.hh"
#include "BasisSpline.hh"
using namespace RavlN;

class BaseStateVectorC
{
	public:
	BaseStateVectorC(){;}
	
	BaseStateVectorC(VectorC const &vect):state(vect.Copy()) {;}
	
	//USEABLE INTERFACE
	virtual VectorC GetStateVector(void) = 0;
	virtual VectorC GetNoiseVector(void) = 0;
	virtual SArray1dC<Point2dC> GetObservationVector(void) = 0;
	virtual BasisSplineC GetBSpline(void) = 0;
	
	protected:
	VectorC state;
};


#endif
