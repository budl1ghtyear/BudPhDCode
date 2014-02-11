#include "BoxPropModel.hh"

ParticleC BoxPropModelC::Propagate(ParticleC &pt)
{
	//This is a zeroth order propagation model
	//x_t = x_t-1 + noise
	//noise = Gaussian random number * adaptive noise window (this is basically this->window*eigenvalue of the state vector)
	//NAIVE PROPAGATION MODEL
	VectorC state = pt.GetState().GetStateVector().Copy();
	//We are now going to provide a new vector N(x_t-1,variance)
	VectorC newstate(state.Size()); newstate.Fill(0.0);
	for(UIntT i = 0; i < state.Size(); i++)
	{
		newstate[i] = GetRandomNumber(state[i],noisewindow);
	}
	pt.GetState().SetStateVector(newstate.Copy());
	//cout<<"Set state is : "<<pt.GetState().GetStateVector()<<endl;
	return pt;	
}

RealT BoxPropModelC::GetRandomNumber(RealT const &mean, RealT const &variance)
{
	//To do this:
	//1) Get a normally distributed random number with 0 mean and unit variance
	RealT gnum = BasePropagationModelC::rnd.Generate();
	RealT result = (gnum*variance*noisewindow)+mean ;
	//cout<<"Random Number Result = "<<result<<"\t with mean = "<<mean<<"\t variance = "<<variance<<"\t window = "<<window<<endl;
	return result;
}

