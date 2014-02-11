#ifndef FOPROPAGATIONMODEL_HH
#define FOPROPAGATIONMODEL_HH
/////////////////////
// File: FOPropagationModelC
// Author: Bud Goswami
// Date: 11th May 2009
// Purpose: Declare a basic first order propagation model, x_t = x_t-1 + NOISE
/////////////////////
#include "BasePropagationModel.hh"
#include "Particle.hh"
#include "BaseStateVector.hh"
#include "Ravl/RandomGauss.hh"
#include "Ravl/Vector2d.hh"

using namespace RavlN;

class FOPropagationModelC : public BasePropagationModelC
{
public:
	FOPropagationModelC(){}
	~ FOPropagationModelC(){}
	FOPropagationModelC(const IntT window);
	//Constructor is for application to a pair of filtered images already (so perform Gaussian Convolution
	ParticleC Propagate(ParticleC &pt);
	
protected:
	static RandomGaussC rnd;
	IntT window;
};

#endif
