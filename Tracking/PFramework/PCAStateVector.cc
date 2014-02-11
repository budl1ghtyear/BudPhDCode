#include "PCAStateVector.hh"
namespace RavlN{
PCAStateVectorBodyC::PCAStateVectorBodyC(DirectoryC const &dir,UIntT const &numcp ):pca(dir,numcp),sim(dir)
{
	VectorC vec(0);
	state = vec.Copy();
}
PCAStateVectorBodyC::PCAStateVectorBodyC(DirectoryC const &dir, VectorC const &pts, UIntT const &numcp):pca(dir,numcp),sim(dir)
{
	state = pts.Copy();
}

PCAStateVectorBodyC::PCAStateVectorBodyC(DirectoryC const &dir, SArray1dC<Point2dC> const &pts, UIntT const &numcp):pca(dir,numcp),sim(dir)
{
	//To go from points in Image space to points in state space requires the following:
	//STEP 1 - projection onto affine space
	Tuple2C<SArray1dC<Point2dC>, VectorC> aff = sim.ProjectToAffineSpace(pts);
	//STEP 2 - projection onto PCA space
	VectorC proj_vect = pca.PerformPCAProjection(aff.Data1());
	//Now concatenate the PCA and affine vectors
	VectorC res = aff.Data2().Append(proj_vect).Copy();
	//~ cout<<"About to append vector to state - "<<res<<endl;
	state = res.Copy();
	//~ cout<<"Copied state is now"<<state<<endl;
}
PCAStateVectorBodyC::PCAStateVectorBodyC(DirectoryC const &dir, Array1dC<Point2dC> const &pts, UIntT const &numcp)
{
	PCAStateVectorBodyC(dir,ArrayToSArray_B(pts),numcp);
}
Array1dC<Point2dC> PCAStateVectorBodyC::GetObservationVector()
{
	//To go from the state vector to the image control points requires the following:
	//Separation into affine and eigen components
	//~ cout<<"State Vector is "<<this->state<<endl;
	VectorC eigen = this->state.After(4);
	cout<<"Eigen part of state - "<<eigen<<endl;
	VectorC affine(this->state[0],this->state[1],this->state[2],this->state[3]);
	SArray1dC<Point2dC> inv_vect = pca.PerformInverseProjection(eigen);
	SArray1dC<Point2dC> resultsarr = sim.ProjectFromAffineSpace(inv_vect,affine);
	return SArrayToArray_B(resultsarr);
}
}
