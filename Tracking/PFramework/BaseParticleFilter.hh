#ifndef BASEPARTICLEFILTERC_HH
#define BASEPARTICLEFILTERC_HH

////////////////////////////////
// Name - BaseParticleFilterC
// Author - Bud Goswami
// Notes : Defines a base particle filter class that has a virtual apply function
////////////////////////////////
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/RefCounter.hh"
#include "Ravl/Stream.hh"
//#include "Ravl/StateVector.hh"
///////////////////////////////
#include "Particle.hh"
#include "BasePropagationModel.hh"
#include "BaseMeasurementModel.hh"
#include "BaseResamplingModel.hh"
#include "BaseStateVector.hh"

using namespace RavlN;

class BaseParticleFilterC
{
public:
	BaseParticleFilterC() : numparticles(0) {}
	virtual ~BaseParticleFilterC() {}
	//: Constructor.
	BaseParticleFilterC(const IntT &n, const BaseStateVectorC &init, const BasePropagationModelC &pmd, const BaseMeasurementModelC &mmd, const BaseResamplingModelC &rmd)
	{
		numparticles = n;
		RealT wt = 1.0/(RealT)(n);
		ParticleC pt(init, wt);
		Array1dC<ParticleC> arr(n);
		for(Array1dIterC<ParticleC> it(arr); it; it++)
		{
			(*it) = pt;
		}
		ptcollection = arr.Copy();
		propmdl = const_cast<BasePropagationModelC*>(&pmd);
		measmdl = const_cast<BaseMeasurementModelC*>(&mmd);
		resampmdl = const_cast<BaseResamplingModelC*>(&rmd);
	}
	BaseParticleFilterC(const IntT &n, const BaseStateVectorC &init, const BasePropagationModelC &pmd, const BaseMeasurementModelC &mmd, const BaseResamplingModelC &rmd, const RealT &defwt)
	{
		numparticles = n;
		ParticleC pt(init, defwt);
		Array1dC<ParticleC> arr(n);
		for(Array1dIterC<ParticleC> it(arr); it; it++)
		{
			(*it) = pt;
		}
		ptcollection = arr.Copy();
		propmdl = const_cast<BasePropagationModelC*>(&pmd);
		measmdl = const_cast<BaseMeasurementModelC*>(&mmd);
		resampmdl = const_cast<BaseResamplingModelC*>(&rmd);
	}
	BaseParticleFilterC(const IntT &n, const BasePropagationModelC &pmd, const BaseMeasurementModelC &mmd, const BaseResamplingModelC &rmd, const Array1dC<ParticleC> &ptc)
	{
		numparticles = n;
		ptcollection = ptc.Copy();
		propmdl = const_cast<BasePropagationModelC*>(&pmd);
		measmdl = const_cast<BaseMeasurementModelC*>(&mmd);
		resampmdl = const_cast<BaseResamplingModelC*>(&rmd);
	}
			
	IntT GetNumParticles(void){return numparticles;}
	Array1dC<ParticleC> GetParticleCollection(void) {return ptcollection;}
	BasePropagationModelC* GetPropagationModel(void) {return propmdl;}
	BaseMeasurementModelC* GetMeasurementModel(void) {return measmdl;}
	BaseResamplingModelC* GetResamplingModel(void) {return resampmdl;}
	BaseStateVectorC GetStepEstimate(void) {return stepestimate;}
	void operator=(const BaseParticleFilterC &other)
	{
		BaseParticleFilterC *oth = const_cast<BaseParticleFilterC*>(&other);
		numparticles = oth->GetNumParticles();
		ptcollection = oth->GetParticleCollection();
		propmdl = oth->GetPropagationModel();
		measmdl = oth->GetMeasurementModel();
		resampmdl = oth->GetResamplingModel();
		stepestimate = oth->GetStepEstimate();
	}
	virtual BaseStateVectorC Apply() = 0;
	
protected:
	
	virtual BaseStateVectorC Estimate()
	{
		VectorC res;
		for(Array1dIterC<ParticleC> it(ptcollection); it; it++)
		{
			VectorC dummy(((*it).GetState().GetStateVector())*((*it).GetWeight()));
			res+= dummy;
		}
		BaseStateVectorC out(res);
		return out;
	}
	IntT numparticles;
	Array1dC<ParticleC> ptcollection;
	BaseStateVectorC stepestimate;
	BasePropagationModelC *propmdl;
	BaseMeasurementModelC *measmdl;
	BaseResamplingModelC *resampmdl;
};

#endif
