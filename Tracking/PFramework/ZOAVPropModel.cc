#include "ZOAVPropModel.hh"

ParticleC ZOAVPropModelC::Propagate(ParticleC &pt)
{
	//This is a zeroth order propagation model
	//x_t = x_t-1 + noise
	//noise = Gaussian random number * adaptive noise window (this is basically this->window*eigenvalue of the state vector)
	//Create noise window
	VectorC noise(pt.GetState().GetNoiseVector(noisewindow));
	for(UIntT i = 0; i < noise.Size(); i++)
	{
		noise[i] *= BasePropagationModelC::rnd.Generate();
	}
	noise += pt.GetState().GetStateVector();
	pt.GetState().SetStateVector(noise);
	return pt;	
}
