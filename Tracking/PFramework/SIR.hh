#ifndef SIR_HH
#define SIR_HH
////////////////////////////////
// Name - SIR
// Author - Bud Goswami
// Notes : Defines a SIR filter
////////////////////////////////
#include "Particle.hh"
#include "BaseParticleFilter.hh"
//#include "Ravl/OS/Date.hh"
#include <ctime>
using namespace std;

class SIR : public BaseParticleFilterC
{
public:
	SIR(){}
	//: Constructor.
	SIR(const IntT &n, const BaseStateVectorC &init, const BasePropagationModelC &pmd, const BaseMeasurementModelC &mmd, const BaseResamplingModelC &rmd)
	: BaseParticleFilterC(n,init,pmd,mmd,rmd) {}
	SIR(const IntT &n, const BaseStateVectorC &init, const BasePropagationModelC &pmd, const BaseMeasurementModelC &mmd, const BaseResamplingModelC &rmd, const RealT &defwt)
	: BaseParticleFilterC(n,init,pmd,mmd,rmd,defwt) {}
	SIR(const IntT &n, const BasePropagationModelC &pmd, const BaseMeasurementModelC &mmd, const BaseResamplingModelC &rmd, const Array1dC<ParticleC> &ptc)
	: BaseParticleFilterC(n,pmd,mmd,rmd,ptc) {}

	virtual BaseStateVectorC Apply();
	virtual BaseStateVectorC Estimate()
	{
		//VectorC res(ptcollection[0].GetState().GetStateVector().Size());
		VectorC res(ptcollection[0].GetStateVector().Size());
		res.Fill(0.0);
		for(Array1dIterC<ParticleC> it(ptcollection); it; it++)
		{
			VectorC dummy(((*it).GetStateVector())*((*it).GetWeight()));
			for(SizeT i = 0; i < res.Size(); i++)
			{
				res[i]+= dummy[i];
			}
		}
		ParticleC res_pt(ptcollection[0]);
		BaseStateVectorC result = res_pt.GetState();
		result.SetStateVector(res);
		//BaseStateVectorC out(res);
		//return res_pt.GetState();
		return result;
	}
};
  
#endif
