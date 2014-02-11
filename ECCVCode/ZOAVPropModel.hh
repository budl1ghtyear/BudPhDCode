#ifndef ZOAVPropModel_HH
#define ZOAVPropModel_HH


#include "BasePropagationModel.hh"

using namespace RavlN;
//Class to implement Zeroth Order Propagation With Adaptive Variance 

class ZOAVPropModelC : public BasePropagationModelC
{
	public:
	ZOAVPropModelC(){}
	ZOAVPropModelC(UIntT const &w=1):noisewindow(w){;}
	ParticleC Propagate(ParticleC &pt);

	protected:
	UIntT noisewindow;
	RealT GetRandomNumber(RealT const &mean, RealT const &variance);
};


#endif
