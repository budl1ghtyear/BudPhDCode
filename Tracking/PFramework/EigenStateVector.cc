#include "EigenStateVector.hh"
namespace RavlN{

EigenStateVectorBodyC::EigenStateVectorBodyC(const Array1dC<Point2dC> &pts, const Affine2dC &headaff, MatrixC &mn, MatrixC &ev, const UIntT &dim)
{
	numeig = dim;
	mean = mn.Copy();
	evect = ev.Copy();
	//Apply Affine Projection To Point Vector
	Array1dC<Point2dC> aff_pts = ProjectAffine(pts,headaff);
	//cout<<"Aff_PTS ="<<aff_pts<<endl;
	head_aff = headaff;
	VectorC pt_vec = ConvertPointArrayToVec(pts);
	//Now do magic stuff with matrices
	evals = ComputeEigenValues(pt_vec, mean, evect).Copy();
	VectorC state_dummy(numeig);
	for(UIntT i = 0; i < numeig; i++)
	{	
		state_dummy[i] = evals[i][0];
	}
	state = state_dummy.Copy();
}

EigenStateVectorBodyC::EigenStateVectorBodyC(const Array1dC<Point2dC> &pts, const Affine2dC &headaff, MatrixC &mn, MatrixC &ev)
{
	numeig = 5; //Standard condition for eigenvalue
	mean = mn;
	evect = ev;
	//Apply Affine Projection To Point Vector
	Array1dC<Point2dC> aff_pts = ProjectAffine(pts,headaff);
	head_aff = headaff;
	VectorC pt_vec = ConvertPointArrayToVec(pts);
	//Now do magic stuff with matrices	
	evals = ComputeEigenValues(pt_vec, mean, evect);
	VectorC state_dummy(numeig);
	for(UIntT i = 0; i < numeig; i++)
	{	
		state_dummy[i] = evals[i][0];
	}
	state = state_dummy.Copy();
}

MatrixC EigenStateVectorBodyC::ComputeEigenValues(const VectorC &pt_vec, const MatrixC &mean, const MatrixC &evect)
{
	MatrixC cp_mat(pt_vec);
	// Lamda = (E)^-1 * (CP - MU);
	MatrixC eig_vals = evect.Inverse()*(cp_mat - mean);
	return eig_vals;
}

Array1dC<Point2dC> EigenStateVectorBodyC::GetObservationVector()
{
	//Perform Reverse Operation
	MatrixC pts = mean + (evect*evals);
	//cout<<"Mean - "<<mean<<"\n Evect*Evals = "<<evect*evals<<endl;
	//cout<<"PTS COMPUTED = "<<pts<<endl;
	SizeT dim = pts.Rows()/2;
	Array1dC<Point2dC> points(dim);
	for(UIntT i =0; i < dim; i++)
	{
		Point2dC pt(pts[i + dim][0], pts[i][0]);
		points[i] = pt.Copy();
	}	
	//Array1dC<Point2dC> answer = ProjectAffine(points,head_aff.Inverse());
	//cout<<"Answer = "<<answer<<endl;
	//return answer;
	return points;
}

Array1dC<Point2dC> EigenStateVectorBodyC::GetPreAffineObservationVector()
{
	//Perform Reverse Operation
	MatrixC pts = mean + (evect*evals);
	SizeT size = pts.Rows()/2;
	Array1dC<Point2dC> points(size);
	for(UIntT i =0; i < size; i++)
	{
		Point2dC pt(pts[i + size][0], pts[i][0]);
		points[i] = pt.Copy();
	}	
	return points;
}
}

