#ifndef BASEMEASUREMENTMODELC_HH
#define BASEMEASUREMENTMODELC_HH

////////////////////////////////
// Name - BaseMeasurementModelC
// Author - Bud Goswami
// Notes : Defines a base measurement model body class that can be used as the base for different child classes
// Follows the example for RCHandleVC derivation
// Implements a reference counted object
////////////////////////////////

#include "Ravl/RefCounter.hh"
#include "Ravl/Stream.hh"
#include "Particle.hh"
using namespace RavlN;


class BaseMeasurementModelC
{
public:
	BaseMeasurementModelC() {}
	virtual ~BaseMeasurementModelC(){}
	//: Constructor.
	virtual RealT Measure(ParticleC &pt) = 0;
	virtual void operator=(const BaseMeasurementModelC &mmdl) {}
};
 
#endif
