#ifndef EIGENSTATEVECTOR_HH
#define EIGENSTATEVECTOR_HH
//! author="Bud Goswami"
//! date="20/8/2008"

#include "BaseStateVector.hh"
#include "Ravl/Matrix3d.hh"
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/Index2d.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/Vector.hh"


namespace RavlN {
  
class EigenStateVectorC : public BaseStateVectorC
{
public:
	EigenStateVector() {} //Default Constructor
	EigenStateVector(const Array1dC<Point2dC> &pts, const Affine2dC &headaff, const MatrixC &mean, const MatrixC &evect);
	
	//Accessor Methods
	EigenValueC GetEigenValues(void) const {return evals;}

protected:
	EigenValueC evals;
};


}

#endif
