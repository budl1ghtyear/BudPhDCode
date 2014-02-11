#include "Particle.hh"

ParticleC::ParticleC()
{
	weight = 0.0;
}
ParticleC::~ParticleC()
{}

ParticleC::ParticleC(const BaseStateVectorC &st,const RealT &wt)
{
	state = st.Copy();
	weight = wt;
}

ParticleC::ParticleC(const BaseStateVectorC &st)
{
	//state = const_cast<BaseStateVectorC*>(&st);
	state = st.Copy();
	weight = 0.0;
}

//Accessor Functions

BaseStateVectorC ParticleC::GetState(void) const
{
	return state;
}

void ParticleC::SetState(const BaseStateVectorC &newst)
{
	//state = const_cast<BaseStateVectorC*>(&newst);
	state = newst;
}
	
RealT ParticleC::GetWeight(void) const
{
	return weight;
}
	
void ParticleC::SetWeight(const RealT &newwt)
{
	weight = newwt;
}

VectorC ParticleC::GetStateVector(void)
{
	return state.GetStateVector();
}
Array1dC<Point2dC> ParticleC::GetObservationVector(void) 
{
	return state.GetObservationVector();
}

