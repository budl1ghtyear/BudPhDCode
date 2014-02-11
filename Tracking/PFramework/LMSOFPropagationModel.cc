#include "LMSOFPropagationModel.hh"

ParticleC LMSOFPropagationModelC::Propagate(ParticleC &pt)
{
	//Get the point array contained in the particle object
	VectorC lip_pts = pt.GetState().GetStateVector().Copy(); 
	//Optical Flow always needs to be used with BSplineStateVectors, user responsibility
	//Add noise to the points
	ImageC<RealT> error_img = of_obj.ErrorVar();
	cout<<"Input State Vector - "<<lip_pts<<endl;
	UIntT size = lip_pts.Size()/2;
	for(UIntT i = 0; i < (size-1); i++)
	{
		//Compute patch Motion
		Index2dC query_ind(lip_pts[i+size],lip_pts[i]);
		cout<<"Index to put in = "<< query_ind << endl;
		ImageRectangleC imrec = this->GetPatchRectangle(query_ind, window);
		Vector2dC mean_mot = ComputePatchMeanMotion(imrec);
		//Error Variance
		/*RealT err_var = 0.0;
		if(error_img.Frame().Contains(query_ind))
			err_var = error_img[query_ind];
		//Error Variation Calculated
		cout<<"Noise Added = " << Sqrt(Abs(err_var))<< endl;
		*/
		RealT gaussiannoise = 5;
		lip_pts[i+size] += (mean_mot.Row() + (rnd.Generate()*gaussiannoise));
		lip_pts[i] += (mean_mot.Col() + (rnd.Generate()*gaussiannoise));
		//lip_pts[i+size] += (mean_mot.Row() + (rnd.Generate()*Sqrt(err_var)));
		//lip_pts[i] += (mean_mot.Col() + (rnd.Generate()*Sqrt(err_var)));
	}
	//ENSURE THAT THE FIRST AND LAST CONTROL POINTS ARE ALWAYS THE SAME
	lip_pts[size-1] = lip_pts[0];
	lip_pts[lip_pts.Size()-1] = lip_pts[size];
	pt.GetState().SetStateVector(lip_pts);	
	cout<<"Output State Vector - "<<lip_pts<<endl;
	//ParticleC out(pt.GetState(),pt.GetWeight());
	return pt;		
}

Vector2dC LMSOFPropagationModelC::ComputePatchMeanMotion(const ImageRectangleC &im)
{
	Vector2dC mean(0.0,0.0);
	//Compute the average optical flow over an image supplied
	if(mot_img.Frame().Contains(im))
	{
		//Compute the mean motion
		ImageC<Vector2dC> patch(mot_img,im);
		for(Array2dIterC<Vector2dC> it(patch); it; it++)
		{
			mean += (*it);
		}
		mean /= patch.Frame().Area();
	}
	else
	{
		cerr<<"The image patch is outside the area of the motion image"<<endl;
	}
	return mean;
}
