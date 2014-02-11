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
using namespace RavlConstN;

static RealT GetGaussianValue(const RealT &sigma, const RealT &xval)
{
	//Generate a 1D gaussian y-value given the value of (x - mu)
	return (1/(sigma*Sqrt(2*pi)))*Exp(-(Sqr(xval)/(2*Sqr(sigma))));
}

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
