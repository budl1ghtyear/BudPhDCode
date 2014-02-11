#ifndef FOPROPAGATIONMODEL_HH
#define FOPROPAGATIONMODEL_HH
/////////////////////
// File: FOPropagationModelC
// Author: Bud Goswami
// Date: 11th May 2009
// Purpose: Declare a basic first order propagation model, x_t = x_t-1 + NOISE
/////////////////////
#include "BasePropagationModel.hh"
#include "Ravl/Vector2d.hh"

using namespace RavlN;

class FOPropagationModelC : public BasePropagationModelC
{
public:
	FOPropagationModelC()
	{
		window = 1;
	}
	FOPropagationModelC(const IntT &w) {window = w;}
	//Constructor is for application to a pair of filtered images already (so perform Gaussian Convolution
	ParticleC Propagate(ParticleC &pt);
	
protected:

	IntT window;
};

#endif
