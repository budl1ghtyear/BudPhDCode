#include "HSOFPropagationModel.hh"


Vector2dC HSOpticalFlowPropagationModelC::ComputeMeanMotion(const ImageC<Vector2dC> &im)
{
	//Process: just mean the whole thing
	Vector2dC sum(0.0,0.0);
	ImageC<Vector2dC> img = im.Copy();
	for(Array2dIterC<Vector2dC> it(img); it; it++)
	{
		sum += (*it).Copy();
	}
	RealT area = (RealT)im.Frame().Area();
	return sum /= area;
}

Vector2dC HSOpticalFlowPropagationModelC::ComputeWeightedMeanMotion(const ImageC<Vector2dC> &im)
{
	Vector2dC sum(0.0,0.0);
	RealT area = 0.0;
	ImageC<Vector2dC> img = im.Copy();
	for(Array2dIterC<Vector2dC> it(img); it; it++)
	{
		if((*it).Modulus() != 0.0)
		{	sum += (*it).Copy();
			area++;
		}
	}
	return sum /= area;	
}

ParticleC HSOpticalFlowPropagationModelC::Propagate(ParticleC &pt)
{
	//Get the point array contained in the particle object
	VectorC lip_pts = pt.GetState().GetStateVector().Copy(); 
	//Optical Flow always needs to be used with PointStateVectors, user responsibility
	//Add noise to the points
	cout<<"Input State Vector - "<<lip_pts<<endl;
	UIntT size = lip_pts.Size()/2;
	for(UIntT i = 0; i < (size-1); i++)
	{
		lip_pts[i+size] += (mean_mot.Row() + (rnd.Generate()*gaussianwindow));
		lip_pts[i] += (mean_mot.Col() + (rnd.Generate()*gaussianwindow));
	}
	//ENSURE THAT THE FIRST AND LAST CONTROL POINTS ARE ALWAYS THE SAME
	lip_pts[size-1] = lip_pts[0];
	lip_pts[lip_pts.Size()-1] = lip_pts[size];
	pt.GetState().SetStateVector(lip_pts);	
	cout<<"Output State Vector - "<<lip_pts<<endl;
	//ParticleC out(pt.GetState(),pt.GetWeight());
	return pt;	
}

