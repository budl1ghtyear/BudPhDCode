#include "EAStateVector.hh"
namespace RavlN{

EAStateVectorBodyC::EAStateVectorBodyC(const Array1dC<Point2dC> &mean, const MatrixC &evect, const Affine2dC &aff, const Array1dC<Point2dC> &ctrl_pts, const Array1dC<RealT> &knot_vec, const UIntT &ord, const UIntT &n_pts):BaseStateVectorBodyC(ConvertPointToVec(ctrl_pts, mean, evect, aff)),BSplineC(ord, knot_vec, ctrl_pts),pca(BaseStateVectorBodyC::ConvertPointArrayToVec(mean), evect)
{
	num_points = n_pts;
}

Array1dC<Point2dC> EAStateVectorBodyC::GetObservationVector()
{
	this->BSplineC::SetControlPoints(this->GetControlPoints().Copy()); //Get the control point vector from the state using the affine transform and eigen transform
	cout<<"INSIDE EASTATEVECTOR CC - Control Points Into BSpline Curve Render are : \n"<< this->BSplineC::GetControlPoints()<<endl;
	return(this->BSplineC::RenderCurve(num_points));
}

Array1dC<Point2dC> EAStateVectorBodyC::GetControlPoints() const
{
	//UIntT aff_dim = state.Size() - 6;
	UIntT aff_dim = 6;
	Affine2dC aff = this->GetAffineTransform();
	//Now just get the EigenValues and get the Control Point Array
	VectorC evals(state.Size() - aff_dim);
	for(UIntT i = 0; i < evals.Size(); i++)
	{
		evals[i] = state[i];
	}
	cout<<"INSIDE EASTATEVECTOR CC - STATE VECTOR IS : \n"<< state<<endl;
	//Now just perform operation
	//VectorC cptsbeforeaffine = pca.Mean() + pca.Projection()*evals; //Potential Issue!
	VectorC cptsbeforeaffine = pca.Mean() + pca.Projection()*evals; //Potential Issue!
	cout<<"INSIDE EASTATEVECTOR CC - MEAN VECTOR IS : \n"<< pca.Mean()<<endl;
	cout<<"INSIDE EASTATEVECTOR CC - CP AFTER EIGEN PROJECTION IS : \n"<< cptsbeforeaffine<<endl;
	//Now perform the operation using the affine transform
	Array1dC<Point2dC> cptsbeforeaffine_arr = BaseStateVectorBodyC::ConvertVecToPointArray(cptsbeforeaffine);
	Array1dC<Point2dC> cpts(cptsbeforeaffine_arr.Size());
	for(Array1dIter2C<Point2dC, Point2dC> it(cptsbeforeaffine_arr,cpts); it; it++)
	{
		Vector2dC vec(it.Data1().Row(),it.Data1().Col());
		//vec = aff.Inverse() * vec;
		//vec = aff * vec;
		it.Data2() = (Point2dC(vec.Row(),vec.Col())).Copy();
	}
	return cpts;
}

Affine2dC EAStateVectorBodyC::GetAffineTransform() const
{
	VectorC aff_val(6),eigvals(state.Size() - 6); // WE KNOW THE WE ARE USING A STRONG AFFINE TRANSFORM AND SO WE USE A 6 PARAMETER VECTOR
	UIntT aff_dim = state.Size() - 6 ;
	Point2dC translation(state[aff_dim],state[aff_dim+1]); cout<<"Translation is now - "<<translation;
	Vector2dC scaling(state[aff_dim+2],state[aff_dim+3]);
	RealT skew = state[aff_dim+4];
	RealT rotation = state[aff_dim+5]; //Now we should be at state.Size) - 1
	Matrix2dC ScalingRotation = Matrix2dC(Cos(rotation),Sin(rotation),-Sin(rotation),Cos(rotation)) * Matrix2dC(scaling[0],skew,skew,scaling[1]) ;
	Affine2dC aff(ScalingRotation, translation);
	cout<<"AFFINE TRANSFORM ACUAL STATE IN CC IS: \n " <<state<<endl;
	cout<<"AFFINE TRANSFORM INSIDE EASTATEVECTOR CC IS: \n " <<aff<<endl;
	return aff;
}
	
}
