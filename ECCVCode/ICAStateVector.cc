#include "ICAStateVector.hh"
namespace RavlN{
ICAStateVectorBodyC::ICAStateVectorBodyC(DirectoryC const &dir):ICA(dir),sim(dir)
{
	VectorC vec(0);
	state = vec.Copy();		
}
ICAStateVectorBodyC::ICAStateVectorBodyC(DirectoryC const &dir, VectorC const &pts):ICA(dir),sim(dir)
{
	state = pts.Copy();
}

ICAStateVectorBodyC::ICAStateVectorBodyC(DirectoryC const &dir, SArray1dC<Point2dC> const &pts):ICA(dir),sim(dir)
{
	//To go from points in Image space to points in state space requires the following:
	//STEP 1 - projection onto affine space
	Tuple2C<SArray1dC<Point2dC>, VectorC> aff = sim.ProjectToAffineSpace(pts);
	//STEP 2 - projection onto ICA space
	VectorC proj_vect = ICA.PerformICAProjection(aff.Data1());
	//Now concatenate the ICA and affine vectors
	VectorC res = aff.Data2().Append(proj_vect).Copy();
	//~ cout<<"About to append vector to state - "<<res<<endl;
	state = res.Copy();
	//~ cout<<"Copied state is now"<<state<<endl;
}
ICAStateVectorBodyC::ICAStateVectorBodyC(DirectoryC const &dir, Array1dC<Point2dC> const &pts)
{
	ICAStateVectorBodyC(dir,ArrayToSArray_B(pts));
}
Array1dC<Point2dC> ICAStateVectorBodyC::GetObservationVector()
{
	//To go from the state vector to the image control points requires the following:
	//Separation into affine and eigen components
	//~ cout<<"State Vector is "<<this->state<<endl;
	VectorC eigen = this->state.From(4);
	cout<<"Eigen part of state - "<<eigen<<endl;
	VectorC affine(this->state[0],this->state[1],this->state[2],this->state[3]);
	SArray1dC<Point2dC> inv_vect = ICA.PerformInverseProjection(eigen);
	SArray1dC<Point2dC> resultsarr = sim.ProjectFromAffineSpace(inv_vect,affine);
	return SArrayToArray_B(resultsarr);
}
}
