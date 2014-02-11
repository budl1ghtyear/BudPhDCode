#include "ZOAVPropModel.hh"

ParticleC ZOAVPropModelC::Propagate(ParticleC &pt)
{
	//This is a zeroth order propagation model
	//x_t = x_t-1 + noise
	//noise = Gaussian random number * adaptive noise window (this is basically this->window*eigenvalue of the state vector)
	//Create noise window
	//~ VectorC noise(pt.GetState().GetNoiseVector(noisewindow));
	//~ VectorC addnoise(noise.Size());addnoise.Fill(0.0);
	//~ VectorC newstate = pt.GetState().GetStateVector().Copy();
	//~ for(UIntT i = 0; i < noise.Size(); i++)
	//~ {
		//~ addnoise[i] = noise[i] * BasePropagationModelC::rnd.Generate() + newstate[i];
		//~ //Now perform verification to ensure that the resulting data is within the valid error bound
		//~ //If the result of noise addition projects the data beyond validity, cap it to the outer bound
		//~ if(addnoise[i] > noise[i])
			//~ addnoise[i] = newstate[i];
	//~ }
	//~ //cout<<"Input Noise Vector = \n"<<noise<<endl;
	//~ pt.GetState().SetStateVector(addnoise.Copy());
	
	//NAIVE PROPAGATION MODEL
	VectorC state = pt.GetState().GetStateVector().Copy();
	VectorC noise = pt.GetState().GetNoiseVector().Copy(); //NOTE THAT THIS SHOULD RETURN THE VARIANCE PARAMETERS!
	//cout<<"NOISE VECTOR IN PROP CLASS = "<<noise<<endl;
	//We are now going to provide a new vector N(x_t-1,variance)
	VectorC newstate(state.Size()); newstate.Fill(0.0);
	for(UIntT i = 0; i < state.Size(); i++)
	{
		if(i < 4 )
		{
			newstate[i] = state[i];
		}
		else
		{ 
			RealT perturbation = GetRandomNumber(state[i],noise[i]);
			newstate[i] = perturbation;
		}
	}
	pt.GetState().SetStateVector(newstate.Copy());
	//cout<<"Set state is : "<<pt.GetState().GetStateVector()<<endl;
	return pt;	
}

RealT ZOAVPropModelC::GetRandomNumber(RealT const &mean, RealT const &variance)
{
	//To do this:
	//1) Get a normally distributed random number with 0 mean and unit variance
	RealT gnum = BasePropagationModelC::rnd.Generate();
	RealT result = (gnum*variance*noisewindow)+mean ;
	//cout<<"Random Number Result = "<<result<<"\t with mean = "<<mean<<"\t variance = "<<variance<<"\t window = "<<window<<endl;
	return result;
}

