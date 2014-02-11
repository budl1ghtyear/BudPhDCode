#include "Particle.hh"

ParticleC::ParticleC()
{
	weight = 0.0;
}
ParticleC::~ParticleC()
{}

ParticleC::ParticleC(const BaseStateVectorC &st,const RealT &wt)
{
	state = st;
	weight = wt;
}

ParticleC::ParticleC(const BaseStateVectorC &st)
{
	state = st;
	weight = 0.0;
}

//Accessor Functions

BaseStateVectorC ParticleC::GetState(void) const
{
	return state;
}

void ParticleC::SetState(const BaseStateVectorC &newst)
{
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
