#include "SystematicResampling.hh"

Array1dC<ParticleC> SystematicResampling::Resample(const Array1dC<ParticleC> &ptc)
{
	Array1dC<ParticleC> ptcol = ptc.Copy();
	IntT numpt = (IntT)ptcol.Size(); //Number of particles
	RealT initwt = (1.0 / (RealT) numpt); //Initial Weight
	Array1dC<RealT> c (ptcol.Size());
	c.Fill(0.0);
	Array1dC<ParticleC> cParticleTmp(ptcol.Size());
	//Calculate the Cumulative Sum of Weights from the Particles
	c[0] = ptcol[0].GetWeight();
	for(IndexC i = 1; i < (IndexC)numpt; i ++)
	{
		c[i] = c[i-1] + ptcol[i].GetWeight();
	}
	//Draw a starting point
	bool inclusive = false; //so random number generator will never return 1
	RealT u1 = (Random1(inclusive))/((RealT)numpt);
	RealT u;
	IndexC uliIndex = 0;
	//Perform resampling
	for(IndexC i = 0; i < ptcol.Size(); i++)
	{
		u = u1 + (i * initwt);
		while((u > c[uliIndex]) && (uliIndex < numpt))
		{
			uliIndex ++;
		}
		cParticleTmp[i] = ptcol[uliIndex];
		cParticleTmp[i].SetWeight(initwt);
	}
	//Copy Back
	for(Array1dIter2C<ParticleC,ParticleC> it(ptcol, cParticleTmp); it; it++)
	{
		it.Data1() = it.Data2();
	}
	return ptcol;
}
RealT SystematicResampling::GetEffiency(const Array1dC<ParticleC> &ptc)
{
	//Perform degenracy algorithm calculation:
	//Neff = 1/Sum of weights squared
	RealT sum = 0.0;
	for(Array1dIterC<ParticleC> it(ptc); it; it++)
	{
		sum += Sqrt((*it).GetWeight());
	}
	efficiency = (1.0 / sum);
	return efficiency;
}
