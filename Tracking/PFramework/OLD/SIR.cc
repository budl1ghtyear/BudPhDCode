#include "SIR.hh"

BaseStateVectorC SIR::Apply()
{
	//RESAMPLE STEP
	ptcollection = resampmdl->Resample(ptcollection);
	//cout<<"STEP 1"<<endl;
	/*
	for(Array1dIterC<ParticleC> it(ptcollection); it; it++)
	{
		cout<<"STATE = "<<(*it).GetState().GetX()<<"\t WEIGHT = "<<(*it).GetWeight()<<endl;
	}*/
	//PROPAGATION AND MEASUREMENT STEP
	for(Array1dIterC<ParticleC> it(ptcollection); it; it++)
	{
		(*it) = propmdl->Propagate((*it));
		//cout<<"Finished Propagation"<<endl;
		(*it).SetWeight(measmdl->Measure((*it)));
		
	}
	//NORMALISE WEIGHTS
	RealT sum = 0.0;
	for(Array1dIterC<ParticleC> it(ptcollection); it; it++)
	{
		sum += (*it).GetWeight();
	}
	for(Array1dIterC<ParticleC> it(ptcollection); it; it++)
	{
		RealT weight = (*it).GetWeight();
		weight /= sum;
		(*it).SetWeight(weight);
	}
	/*
	cout<<"STEP 2"<<endl;
	for(Array1dIterC<ParticleC> it(ptcollection); it; it++)
	{
		cout<<"STATE = "<<(*it).GetState().GetX()<<"\t WEIGHT = "<<(*it).GetWeight()<<endl;
	}*/	
	//ESTIMATE
	stepestimate = Estimate();
	//cout<<"FINAL ESTIMATE\t "<<stepestimate.GetX()<<endl;
	return stepestimate;
}
