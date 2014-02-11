#ifndef BASEPROPAGATIONMODELC_HH
#define BASEPROPAGATIONMODELC_HH

////////////////////////////////
// Name - BasePropagationModelC
// Author - Bud Goswami
// Notes : Defines a base propagation model body class that can be used as the base for different child classes
// Follows the example for RCHandleVC derivation
// Implements a reference counted object
////////////////////////////////

#include "Ravl/RefCounter.hh"
#include "Ravl/Stream.hh"
#include "Particle.hh"
using namespace RavlN;

class BasePropagationModelC 
{
public:
	BasePropagationModelC() {}
	virtual ~BasePropagationModelC(){}
	//: Constructor.
	virtual ParticleC Propagate(ParticleC &pt) = 0;
	virtual void operator=(const BasePropagationModelC &mmdl)
	{}
};

#endif
