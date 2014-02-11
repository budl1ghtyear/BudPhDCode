#include "LLMeasurementModel.hh"

LLMeasurementModelC::LLMeasurementModelC(const ImageC<RealT> &jim, const MeanCovarianceC &lp, const MeanCovarianceC &sk, const UIntT &win)
{
	vec_img = jim;
	lip = lip;
	skin = sk;
	window = win;	
}
LLMeasurementModelC::LLMeasurementModelC(const ImageC<RealT> &jim, const MeanCovarianceC &lp, const MeanCovarianceC &sk)
{
	vec_img = jim;
	lip = lip;
	skin = sk;
	window = 10;	//default search window of 10 pixels
}
RealT LLMeasurementModelC::Measure(ParticleC &pt)
{
	//Array1dC<Point2dC> ind = pt.GetState().GetObservationVector();
	Array1dC<Point2dC> ind = pt.GetObservationVector();
	//cout<<"Observation Vector "<<ind;
	RealT weight = ComputeWeight(ind,(pt.GetState()).GetBSpline());
	return weight;	
}

//RealT LLMeasurementModelC::GetGaussianValue(const RealT &sigma, const RealT &xval)
//{
//	//Generate a 1D gaussian y-value given the value of (x - mu)
//	return (1/(sigma*Sqrt(2*pi)))*Exp(-(Sqr(xval)/(2*Sqr(sigma))));
//}

RealT LLMeasurementModelC::ComputeWeight(const Array1dC<Point2dC> &ar)
{
	//Using Points, Get B-Spline Control Points
	UIntT order = 3; UIntT ncp = 11;
	BSplineC bspl(order,ncp, BSplineC::UOPEN);
	bspl.SetControlPoints(ar);
	return doComputeWeight(ar,bspl);
}
RealT LLMeasurementModelC::ComputeWeight(const Array1dC<Point2dC> &ar, const BSplineC &bspl)
{
	return doComputeWeight(ar,bspl);
}

RealT LLMeasurementModelC::doComputeWeight(const Array1dC<Point2dC> &ar, const BSplineC &bspl)
{
	//Now Measure the distance
	UIntT curve_res = 20;
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
		//for(RealT i = -5; i < search; i++)
		for(RealT i = search; i > 0; i--)
		{
			//Point2dC cur = normal.Point(i); Point2dC next = normal.Point(i+1);
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
