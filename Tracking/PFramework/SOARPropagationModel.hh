#ifndef SOARPROPAGATIONMODEL_HH
#define SOARPROPAGATIONMODEL_HH
/////////////////////
// File: FOPropagationModelC
// Author: Bud Goswami
// Date: 11th May 2009
// Purpose: Declare a basic first order propagation model, x_t = x_t-1 + NOISE
/////////////////////
#include "BasePropagationModel.hh"
#include "Ravl/Matrix.hh"

using namespace RavlN;

class SOARPropagationModelC : public BasePropagationModelC
{
public:
	SOARPropagationModelC(){}
	SOARPropagationModelC(const IntT &w) {window = w;}
	SOARPropagationModelC(const MatrixC &a_one, const MatrixC &a_two, const MatrixC &dist, const VectorC &prev_state, const IntT &wind)
	{
		a1 = a_one;
		a2 = a_two;
		d = dist;
		window = wind;	
		prev = prev_state;
	}
	SOARPropagationModelC(const MatrixC &a_one, const MatrixC &a_two, const MatrixC &dist, const VectorC &prev_state)
	{
		a1 = a_one;
		a2 = a_two;
		d = dist;
		prev = prev_state;
		window = 1;	//default value of the gaussian noise window. This is set to 1 so we don't get any scaling in the noise
	}
	//Constructor is for application to a pair of filtered images already (so perform Gaussian Convolution
	ParticleC Propagate(ParticleC &pt);
protected:
	
	IntT window;
	MatrixC a1;
	MatrixC a2;
	MatrixC d;
	VectorC prev;
};

#endif
