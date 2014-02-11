#ifndef AFFINEPOINTSTATEVECTOR_HH
#define AFFINEPOINTSTATEVECTOR_HH
//! author="Bud Goswami"
//! date="17/05/2009"

#include "PointStateVector.hh"
//This class is to deal with affine normalised euclidean points

namespace RavlN {
  
class AffinePointStateVectorC : public PointStateVectorC
{
public:
	AffinePointStateVectorC(){}
	//CONSTRUCTORS EXPECT THE POINTS TO ALREADY BE AFFINE NORMALISED
	AffinePointStateVectorC(Array1dC<Point2dC> &ar, const Affine2dC &af):PointStateVectorC(ar){head_aff = af;}
	AffinePointStateVectorC(Array1dC<Index2dC> &ar, const Affine2dC &af):PointStateVectorC(ar){head_aff = af;}
	
	AffinePointStateVectorC(const AffinePointStateVectorC &pt)
	{
		state = pt.state.Copy();
	}
	Array1dC<Point2dC> GetPreAffineObservationVector()
	{
		return PointStateVector::GetObservationVector();
	}
	
	Array1dC<Point2dC> GetObservationVector()
	{
		Array1dC<Point2dC> aff_ar = ProjectAffine(PointStateVector::GetObservationVector(),head_aff);
		return aff_ar;
	}
	//Accessor Methods:
	Affine2dC GetAffine(void) const {return head_aff;}
	void SetAffine(const Affine2dC &af) {head_aff = af;}
	
protected:
	Affine2dC head_aff;
};


}

#endif
