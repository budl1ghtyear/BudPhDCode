#ifndef PARTICLEC_HH
#define PARTICLEC_HH

/////////////////////
// File: ParticleC.hh
// Author: Bud Goswami
// Date: 18th August 2008
// Purpose: Declare the particle class for use in particle filters
/////////////////////

#include "BaseStateVector.hh"
using namespace RavlN;

class ParticleC
{
public:
	ParticleC();
	~ParticleC();
	ParticleC(const BaseStateVectorC &st,const RealT &wt);
	ParticleC(const BaseStateVectorC &st);
	//Accessor Functions
	BaseStateVectorC GetState(void) const;
	void SetState(const BaseStateVectorC &newst);
	void SetStateVector(const VectorC &vec)
	{
		state.SetStateVector(vec);
	}
	VectorC GetStateVector(void);
	Array1dC<Point2dC> GetObservationVector(void);
	
	RealT GetWeight(void) const;
	void SetWeight(const RealT &newwt);
	void operator=(const ParticleC &oth)
	{
		state = oth.state.Copy();
		weight = oth.GetWeight();
	}
	//~ //Copy Constructor
	ParticleC(const ParticleC &pt)
	{
		state= pt.GetState().Copy();
		//pt.GetState().Copy();
		weight = pt.GetWeight();
	}
	ParticleC Copy() const
	{
		return ParticleC(state,weight);
	}
protected:
	//Data Members
	BaseStateVectorC state;
	RealT weight;
};

#endif
