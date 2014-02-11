#ifndef BOXPropModel_HH
#define BOXPropModel_HH


#include "BasePropagationModel.hh"

using namespace RavlN;
//Class to implement Zeroth Order Propagation With Adaptive Variance 

class BoxPropModelC : public BasePropagationModelC
{
	public:
	BoxPropModelC(){}
	BoxPropModelC(UIntT const &w=1):noisewindow(w){;}
	ParticleC Propagate(ParticleC &pt);

	protected:
	UIntT noisewindow;
	RealT GetRandomNumber(RealT const &mean, RealT const &variance);
};

#endif
