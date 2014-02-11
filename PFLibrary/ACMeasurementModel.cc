#include "ACMeasurementModel.hh"

RealT ACMeasurementModelC::Measure(ParticleC &pt)
{
	SArray1dC<Point2dC> ind = ArrayToSArray_B(pt.GetObservationVector());
	RealT weight = lamda * ComputeInternalEnergy(ind) + ComputeExternalEnergy(ind,(pt.GetState()).GetBSpline()); //formula for the energy term which is the measurement
	return weight;	
}

RealT ACMeasurementModelC::ComputeInternalEnergy(const SArray1dC<Point2dC> &cur) const 
{
	RealT answer = 0.0;
	//This is the internal energy term. It is the same equation as described in Section 5.1 of the paper
	//For ease of use, since this equation requires wrapping around, we are going to use the DListC class
	DListC<Point2dC>  cur_data;
	for(SArray1dIterC<Point2dC> it(cur); it; it++)
	{
		cur_data.Append((*it).Copy());
	} 
	//Now iterate through the DList and solve the equations
	for(UIntT i = 0; i < cur_data.Size(); i++)
	{
		if(i == cur_data.Size() - 1)
		{
			answer += Point2dC(cur_data.Nth(i-1) - cur_data.Nth(i)).EuclideanDistance(Point2dC(cur_data.Nth(i) - cur_data.Nth(0)));
		}
		else
		{
			answer += Point2dC(cur_data.Nth(i-1) - cur_data.Nth(i)).EuclideanDistance(Point2dC(cur_data.Nth(i) - cur_data.Nth(i+1)));
		}
	}
	cout<<"Computed Internal Energy"<<endl;
	return answer;
}

RealT ACMeasurementModelC::ComputeExternalEnergy(const SArray1dC<Point2dC> &cur, const BasisSplineC &bspl)  
{
	//RealT answer = lamda * ComputeContourEnergy(cur,bspl) + (1 - lamda) * ComputeMotionCompensationTerm(cur);
	RealT answer = ComputeContourEnergy(cur,bspl);
	return answer;
}
RealT ACMeasurementModelC::ComputeContourEnergy(const SArray1dC<Point2dC> &cur, const BasisSplineC &bspl) 
{
	//We have been provided an edge image of some description
	//Now Measure the distance
	SArray1dC<LinePP2dC> normals = bspl.GetCurveNormals();
	if(normals.Size() == 0)
	{
		cerr<<"BSpline has not been initialised correctly"<<endl;
		exit(1);
	}
	//Now we want to go around the curve measuring the distance from the curve to the nearest edge
	RealT j_wt = 0.0;IntT search = window;
	for(UIntT pval = 0; pval < normals.Size(); pval ++)		
	{
		LinePP2dC normal = normals[pval]; //Get B-Spline Normal Line
		Point2dC start = normal.Point(0.0);
		//Now just search from the start to the end for the binary condition to end
		//Change condition to start looking outwards in and search for the maximum value which corresponds to the strongest edge
		RealT index = 0.0; 
		RealT maxval = 0.0;
		for(RealT i = search; i > 0; i--)
		{
			Point2dC cur = normal.Point(i); 
			Index2dC current(cur.Row(),cur.Col());
			if(vec_img.Contains(current))
			{
				if(vec_img[current] > maxval)
				{
					//to get here: vec_img must contain the generated index and the value at that point must be larger than the previous value of maxval
					index = search;
				}
			}
		}
		j_wt += GetGaussianValue(3,index);
	}
	cout<<"Computed Contour Energy"<<endl;
	return j_wt;
}
	
RealT ACMeasurementModelC::ComputeMotionCompensationTerm(const SArray1dC<Point2dC> &cur) 
{
	//Compute Motion Compensation Term as described in Section 4
	IntT mcwindow = 8;
	SArray1dC<Vector2dC> displacementVectors(cur.Size());
	SArray1dC<RealT> displacements(cur.Size());
	RealT answer = 0.;
	for (SArray1dIter3C<Point2dC,Vector2dC,RealT> it(cur,displacementVectors,displacements);it;it++)
	{
		it.Data2() = Vector2dC(0.0, 0.0);
		
		ImageRectangleC pointPatchCheck( it.Data1().Row()-mcwindow, it.Data1().Row()+mcwindow, it.Data1().Col()-mcwindow, it.Data1().Col()+mcwindow);
		ImageRectangleC toCheck(oflow_img.Frame());
		ImageRectangleC pointPatch = RangeCheck(pointPatchCheck,toCheck);
		cout<<"PointPatchCheck - "<<pointPatchCheck<<"\t PointPatch - "<<pointPatch<<"\t toCheck - "<<toCheck<<endl;
		for (Array2dIterC<Vector2dC> it2(oflow_img,pointPatch);it2;it2++)
		{
			it.Data2() += (*it2);
		}
		it.Data2() /= ( (2 * window+ 1) * (2 * window + 1) );
		it.Data3() = it.Data2().Norm();
		answer += it.Data3();
	}
	cout<<"Computed Motion Compensation Energy"<<endl;
	return answer;
}
