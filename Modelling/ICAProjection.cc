#include "ICAProjection.hh"

ICAProjectionC::ICAProjectionC(DirectoryC const &dir)
{
	FilenameC mix = dir+"Amatrix.tsv";
	FilenameC umix = dir+"Wmatrix.tsv";
	if((!LoadFileFromStream(mix, this->mixing))||(!LoadFileFromStream(umix,this->unmixing)))
	{
		cerr<<"Could not load the ICA files"<<endl;
		exit(1);
	}	
	cout<<"INST the class"<<endl;
}
	
VectorC ICAProjectionC::PerformICAProjection(SArray1dC<Point2dC> const &points)
{
	VectorC vect = PointsToVector(points);
	if(vect.Size() != unmixing.Cols())
	{
		cerr<<"Matrix Inner Dimensions MUST agree for multiplication"<<endl;
		exit(1);
	}
	VectorC res = this->unmixing * vect;
	return res;
}
VectorC ICAProjectionC::PerformICAProjection(Array1dC<Point2dC> const &points)
{
	return VectorC(PerformICAProjection(ArrayToSArray_B(points)).Copy());
}
SArray1dC<Point2dC> ICAProjectionC::PerformInverseProjection(VectorC const &vect)
{
	if(vect.Size() != mixing.Cols())
	{
		cerr<<"Matrix Inner Dimensions MUST agree for multiplication"<<endl;
		exit(1);
	}	
	VectorC pts = this->mixing * vect;
	return SArray1dC<Point2dC>(VectorToPoints(pts).Copy());
}
