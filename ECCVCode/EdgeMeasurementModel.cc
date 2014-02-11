#include "EdgeMeasurementModel.hh"

EdgeMeasurementModelC::EdgeMeasurementModelC(ImageC<RealRGBValueC> const &img, ImageRectangleC const &imrec, IntT const &win):im(img.Copy()),myimgrect(imrec),searchwindow(win)
{
	edge = GetSusanImage(ConvertRGBToRealImage(im)).Copy();
	//if(!Save("@X: SUSAN IMAGE",edge)) cerr<<"Could not show SUSAN image"<<endl;
}
RealT EdgeMeasurementModelC::Measure(ParticleC &pt)
{
	//Procedure for measurement is:
	//Get B-Spline Points and Normals from the state vector
	//For each normal, perform a search for the maximum edge in the normal direction within the search window
	//the weight is the Gaussian weighted summation of the distance to the maximum edge for all rendered points
	
	BasisSplineC spl(ArrayToSArray_B(pt.GetState().GetObservationVector()));
	SArray1dC<LinePP2dC> norms = spl.GetCurveNormals();
	return this->MeasureQuality(norms);
	//~ RealT weight = 0.0;
	//~ //iterate through all the rendered points
	//~ for(SArray1dIterC<LinePP2dC> it(norms); it; it++)
	//~ {
		//~ //Look for max edge within searchwindow number of pixels along the normal
		//~ //cout<<"\n New edge measurement : "<<endl;
		//~ RealT maxval = 100.0; RealT dist=0.0;
		//~ for(IntT i = -searchwindow; i < searchwindow; i++)
		//~ {
			//~ Point2dC pt = (*it).Point((RealT)i);
			//~ Index2dC cur(pt.Row(),pt.Col());
			//~ if((edge.Contains(cur))&&(myimgrect.Contains(cur))) //if the edge is not inside the image, its a bad particle and we have a situation where the pixel isn't even considered
			//~ {
				//~ if(edge[cur] > maxval)
				//~ {
					//~ maxval = edge[cur];
					//~ dist = i;
				//~ }
				//~ //cout<<edge[cur]<<"\t";
			//~ }
			//~ //cout<<cur<<"\t";
		//~ }
		//~ //cout<<"\n";
		//~ //cout<<"Maximum distance is : "<<dist<<endl;
		//~ //At this point, the value of dist should be the index from the current rendered point at which the maximum edge value if obtained
		//~ weight += GetGaussianValue(3,dist);
		//~ //cout<<"Weight is now: "<<weight<<endl;
	//~ }
	//~ //Assign the value of the weight to the Particle in question
	//~ return weight;
}
RealT EdgeMeasurementModelC::MeasureQuality(SArray1dC<LinePP2dC> &norms)
{	
	RealT weight = 0.0;
	//iterate through all the rendered points
	for(SArray1dIterC<LinePP2dC> it(norms); it; it++)
	{
		//Look for max edge within searchwindow number of pixels along the normal
		//cout<<"\n New edge measurement : "<<endl;
		RealT maxval = 100.0; RealT dist=0.0;
		for(IntT i = -searchwindow; i < searchwindow; i++)
		{
			Point2dC pt = (*it).Point((RealT)i);
			Index2dC cur(pt.Row(),pt.Col());
			if((edge.Contains(cur))&&(myimgrect.Contains(cur))) //if the edge is not inside the image, its a bad particle and we have a situation where the pixel isn't even considered
			{
				if(edge[cur] > maxval)
				{
					maxval = edge[cur];
					dist = i;
				}
				//cout<<edge[cur]<<"\t";
			}
			//cout<<cur<<"\t";
		}
		//cout<<"\n";
		//cout<<"Maximum distance is : "<<dist<<endl;
		//At this point, the value of dist should be the index from the current rendered point at which the maximum edge value if obtained
		weight += GetGaussianValue(3,dist);
		//cout<<"Weight is now: "<<weight<<endl;
	}
	//Assign the value of the weight to the Particle in question
	return weight;
}
