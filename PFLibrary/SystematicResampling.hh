#ifndef SYSTEMATICRESAMPLINGMODELC_HH
#define SYSTEMATICRESAMPLINGMODELC_HH

////////////////////////////////
// Name - SystematicResamplingModelC
// Author - Bud Goswami
// Notes : Defines a systematic resampling model body class that can be used as the base for different child classes
// Implements Ristic code for systematic resampling
////////////////////////////////

#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/Array1dIter2.hh"
#include "Ravl/Random.hh"
#include "Ravl/RefCounter.hh"
#include "Ravl/Stream.hh"
#include "BaseResamplingModel.hh"
using namespace RavlN;


class SystematicResamplingC : public BaseResamplingModelC 
{
public:
	SystematicResamplingC() {};
	~SystematicResamplingC() {};

	RealT GetEffiency(const Array1dC<ParticleC> &ptc);
	RealT GetEffMeasure(void) const{return efficiency;}
	void SetEfficiency(const RealT &eff) {efficiency = eff;} 
	Array1dC<ParticleC> Resample(const Array1dC<ParticleC> &ptc);
	virtual void operator=(const SystematicResamplingC &mmdl) 
	{
		efficiency = mmdl.GetEffMeasure();
	}
private:
	RealT efficiency;
};
  
#endif
