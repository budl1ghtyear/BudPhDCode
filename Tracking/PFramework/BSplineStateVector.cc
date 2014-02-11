#include "BSplineStateVector.hh"
namespace RavlN {
BSplineStateVectorBodyC::BSplineStateVectorBodyC(const Array1dC<Point2dC> &control_pts, const Array1dC<RealT> &knot_vec, const UIntT &num_pts, const UIntT &ord):BaseStateVectorBodyC(ConvertPointToVec(control_pts)),BSplineC(ord, knot_vec,control_pts)//bSpline(ord, knot_vec,control_pts)
{
	//ctrl_points = control_pts.Copy();	
	//k_vector = knot_vec.Copy();	
	num_points = num_pts;
	//order = ord;	
}
Array1dC<Point2dC> BSplineStateVectorBodyC::GetObservationVector()
{
	return(BSplineC::RenderCurve(num_points));
	//return(bSpline.RenderCurve(num_points));
	
}
}
