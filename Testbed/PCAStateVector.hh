#ifndef PCAStateVector_HH
#define PCAStateVector_HH

#include "BaseStateVector.hh"
#include "PCAProjection.hh"
#include "SimilarityProjection.hh"
using namespace RavlN;

class PCAStateVectorC : public BaseStateVectorC
{
	PCAStateVectorC(){;}
	PCAStateVectorC(VectorC const &pcavect,DirectoryC const &pcadir, UIntT const &numeigen);
	PCAStateVectorC(SArray1dC<Point2dC> const &ctrlpts,DirectoryC const &pcadir, UIntT const &numeigen);
	VectorC GetStateVector(void);
	VectorC GetNoiseVector(void);
	SArray1dC<Point2dC> GetObservationVector(void);
	BasisSplineC GetBSpline(void);	
	protected:
	PCAProjectionC pca;
	SimilarityProjectionC sim;
}

#endif
