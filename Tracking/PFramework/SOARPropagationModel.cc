#include "SOARPropagationModel.hh"

ParticleC SOARPropagationModelC::Propagate(ParticleC &pt)
{
	//Assuming the SOAR Parameter Estimation has already been completed, calculate using:
	//S_n = A2*S_n-2 + A1*S_n-1 + D(perturbation vector)
	//VectorC cur_st = pt.GetState().GetStateVector();
	VectorC cur_st = pt.GetStateVector();
	MatrixC theta_prev(prev);
	MatrixC theta_cur(cur_st);
	MatrixC noise(d);
	for(UIntT i = 0; i < noise.Rows(); i ++)
	{
		noise[i][0] += (BasePropagationModelC::rnd.Generate()*window);
	}
	MatrixC theta = a2*theta_prev + a1*theta_cur + noise;
	VectorC res(theta.Rows());
	for(UIntT i = 0; i < res.Size(); i++)
	{
		res[i] = theta[i][0];
	}
	cout<<res<<endl;
	pt.GetState().SetStateVector(res);
	return pt;
}
