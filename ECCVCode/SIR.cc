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
	//cout<<"INIT STATE = "<<ptcollection[0].GetState().GetStateVector()<<"\t WEIGHT = "<<ptcollection[0].GetWeight()<<endl;
	//PROPAGATION AND MEASUREMENT STEP
	//Can be implemented in same loop, just doing this for speed testing purposes
	time_t systime;
	systime = time(NULL);
	//cout<<"Start time"<<ctime(&systime)<<"\t";
	for(Array1dIterC<ParticleC> it(ptcollection); it; it++)
	{
		(*it) = propmdl->Propagate((*it));
		//~ cout<<"Propagated particle - "<<(*it).GetState().GetStateVector()<<endl;
	}
	systime = time(NULL);
	//cout<<"Finished Propagation \t Time - "<<ctime(&systime)<<endl;
	for(Array1dIterC<ParticleC> it(ptcollection); it; it++)
	{
		//cout<<"Propagated observation - "<<(*it).GetState().GetObservationVector()<<endl;
		(*it).SetWeight(measmdl->Measure((*it)));
		//~ cout<<"Obtained weight is: "<<it->GetWeight()<<endl;
	}
	systime = time(NULL);
	//~ cout<<"Finished Measurement \t Time - "<<ctime(&systime)<<endl;	
	
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
	//~ cout<<"About to perform estimation"<<endl;
	//ESTIMATE
	stepestimate = Estimate();
	//cout<<"FINAL ESTIMATE\t "<<stepestimate.GetStateVector()<<endl;
	return stepestimate;
}

BaseStateVectorC SIR::ApplyGroup()
{
	//RESAMPLE STEP
	ptcollection = resampmdl->Resample(ptcollection);
	//cout<<"INIT STATE = "<<ptcollection[0].GetState().GetStateVector()<<"\t WEIGHT = "<<ptcollection[0].GetWeight()<<endl;
	//PROPAGATION STEP
	//Can be implemented in same loop, just doing this for speed testing purposes
	for(Array1dIterC<ParticleC> it(ptcollection); it; it++)
	{
		(*it) = propmdl->Propagate((*it));
		//cout<<"Propagated particle - "<<(*it).GetState().GetStateVector()<<endl;
	}
	//Group particle measurement
	SampleC<SArray1dC<Point2dC> > curvepoints;
	SampleC<SArray1dC<LinePP2dC> > curvenormals;
	GroupBSpline(ptcollection,curvepoints,curvenormals);
	for(UIntT i = 0; i < ptcollection.Size(); i++)
	{
		ptcollection[i].SetWeight(measmdl->MeasureQuality(curvenormals[i]));
	}
	//NORMALISE WEIGHTS
	RealT sum = 0.0;
	for(Array1dIterC<ParticleC> it(ptcollection); it; it++)
	{
		sum += (*it).GetWeight();
	}
	//~ cout<<"SIR: Sum is "<<sum<<endl;
	for(Array1dIterC<ParticleC> it(ptcollection); it; it++)
	{
		RealT weight = (*it).GetWeight();
		weight /= sum;
		(*it).SetWeight(weight);
	}
	//ESTIMATE
	stepestimate = Estimate();
	cout<<"FINAL ESTIMATE\t "<<stepestimate.GetStateVector()<<endl;
	return stepestimate;
}
