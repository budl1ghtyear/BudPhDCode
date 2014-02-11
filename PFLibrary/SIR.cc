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
	//Can be implemented in same loop, just doing this for speed testing purposes
	time_t systime;
	systime = time(NULL);
	//cout<<"Start time"<<ctime(&systime)<<"\t";
	for(Array1dIterC<ParticleC> it(ptcollection); it; it++)
	{
		//cout<<"Propagated particle - "<<(*it).GetState().GetStateVector()<<endl;
		(*it) = propmdl->Propagate((*it));
	}
	systime = time(NULL);
	//cout<<"Finished Propagation \t Time - "<<ctime(&systime)<<endl;
	for(Array1dIterC<ParticleC> it(ptcollection); it; it++)
	{
		//cout<<"Propagated observation - "<<(*it).GetState().GetObservationVector()<<endl;
		(*it).SetWeight(measmdl->Measure((*it)));
	}
	systime = time(NULL);
	//cout<<"Finished Measurement \t Time - "<<ctime(&systime)<<endl;	
	
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
	//ESTIMATE
	stepestimate = Estimate();
	//cout<<"FINAL ESTIMATE\t "<<stepestimate.GetX()<<endl;
	return stepestimate;
}
