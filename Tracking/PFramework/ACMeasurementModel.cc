#include "ACMeasurementModel.hh"

ACMeasurementModel::ACMeasurementModel(const ImageC<RealT> &like_vec_img, const ImageC<Vector2dC> &img, const Array1dC<Point2dC> &prev, const RealT &g = 0.5, const RealT &l = 0.5, const RealT &b = 1, const UIntT &w = 8)
{
	vec_img = like_vec_img.Copy();
	oflow_img = img.Copy();
	prev_snaxels = prev.Copy();
	gamma = g;
	lamda = l;
	beta = b;
	window = w;
}

RealT ACMeasurementModel::Measure(ParticleC &pt)
{
	Array1dC<Point2dC> ind = pt.GetObservationVector();
	RealT weight = lamda * ComputeInternalEnergy(ind) + ComputeExternalEnergy(ind,(pt.GetState()).GetBSpline()); //formula for the energy term which is the measurement
	return weight;	
}

RealT ACMeasurementModel::ComputeInternalEnergy(const Array1dC<Point2dC> &cur) const 
{
	RealT answer = 0.0;
	//This is the internal energy term. It is the same equation as described in Section 5.1 of the paper
	//For ease of use, since this equation requires wrapping around, we are going to use the DListC class
	DListC<Point2dC>  cur_data;
	for(Array1dC<Point2dC> it(cur); it; it++)
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
	return answer;
}

RealT ACMeasurementModel::ComputeExternalEnergy(const Array1dC<Point2dC> &cur, const BSplineC &bspl) const 
{
	RealT answer = lamda * ComputeContourEnergy(cur,bspl) + (1 - lamda) * ComputeMotionCompensationTerm(cur);
	return answer;
}
RealT ACMeasurementModel::ComputeContourEnergy(const Array1dC<Point2dC> &cur, const BSplineC &bspl) const
{
	//We have been provided an edge image of some description
	//Now Measure the distance
	UIntT curve_res = cur.Size();
	RealT incr = 1.0/(RealT)curve_res;	
	RealT j_wt = 0.0;IntT search = window;
	for(RealT pval = 0; pval < curve_res; pval ++)		
	{
		RealT par = (RealT)pval*incr;
		LinePP2dC normal = bspl.GetNormalLine(par); //Get B-Spline Normal Line
		Point2dC start = normal.Point(0.0);
		//Now just search from the start to the end for the binary condition to end
		//Change condition to start looking outwards in
		RealT index = 0.0;
		for(RealT i = search; i > 0; i--)
		{
			Point2dC cur = normal.Point(i); Point2dC next = normal.Point(i-1);
			if((vec_img.Frame().Contains(next)) && (vec_img.Frame().Contains(cur)))
			{
				if((vec_img[cur] <= 0)&&(vec_img[next] > 0))
				{
					index = i;
					break;
				}
			}
		}
		j_wt += GetGaussianValue(3,index);
	}
	return j_wt;
}
	
RealT ACMeasurementModel::ComputeMotionCompensationTerm(const Array1dC<Point2dC> &cur) const
{
	//Compute Motion Compensation Term as described in Section 4
	Array1dC<Vector2dC> displacementVectors(cur.Size());
	Array1dC<RealT> displacements(cur.Size());
	RealT answer = 0.;
	for (Array1dIter3C<Point2dC,Vector2dC,RealT> it(cur,displacementVectors,displacements);it;it++)
	{
		it.Data2() = Vector2dC(0.0, 0.0);
		ImageRectangleC pointPatch( it.Data1().Row()-window, it.Data1().Row()+window, it.Data1().Col()-window, it.Data1().Col()+window);
		for (Array2dIterC<Vector2dC> it2(oflow_img,pointPatch);it2;it2++)
		{
			it.Data2() += (*it2);
		}
		it.Data2() /= ( (2 * winSize + 1) * (2 * winSize + 1) );
		it.Data3() = it.Data2().Norm();
		answer += it.Data3();
	}
	return answer;
}
