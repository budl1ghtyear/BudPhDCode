#include "FOPropagationModel.hh"

ParticleC FOPropagationModelC::Propagate(ParticleC &pt)
{
	//for the particle, first obtain the state
	VectorC prev_state = pt.GetState().GetStateVector();
	//Now simply add noise to the state vector
	for(UIntT i = 0; i < prev_state.Size(); i++)
	{
		prev_state[i] += (BasePropagationModelC::rnd.Generate()*window); //current state = previous value perturbed by some window 
	}
	cout<<prev_state<<endl;
	pt.GetState().SetStateVector(prev_state);
	return pt;
}
