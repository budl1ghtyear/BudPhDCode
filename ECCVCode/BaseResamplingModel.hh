#ifndef BASERESAMPLINGMODELC_HH
#define BASERESAMPLINGMODELC_HH

////////////////////////////////
// Name - BaseResamplingModelC
// Author - Bud Goswami
// Notes : Defines a base resampling model body class that can be used as the base for different child classes
// Follows the example for RCHandleVC derivation
// Implements a reference counted object
////////////////////////////////

#include "Ravl/Array1d.hh"
#include "Ravl/RefCounter.hh"
#include "Ravl/Stream.hh"
#include "Particle.hh"
using namespace RavlN;
//! userlevel=Develop
//: Example body base class.

class BaseResamplingModelC 
{
public:
	BaseResamplingModelC() {}
	virtual ~BaseResamplingModelC(){}
	//: Constructor.
	virtual Array1dC<ParticleC> Resample(const Array1dC<ParticleC> &ptc) = 0;
	virtual void operator=(const BaseResamplingModelC &mmdl) {}
	
};
  
#endif
