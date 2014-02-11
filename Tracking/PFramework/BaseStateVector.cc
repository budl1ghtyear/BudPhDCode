#include "BaseStateVector.hh"

namespace RavlN
{
Array1dC<Point2dC> BaseStateVectorBodyC::ConvertVecToPointArray(const VectorC &v) const
{
	UIntT size = v.Size()/2;
	Array1dC<Point2dC> out(size);
	for(UIntT i = 0; i < size; i++)
	{	
		Point2dC pt(v[i + size], v[i]);
		out[i] = pt.Copy();
	}
	return out;
}

VectorC BaseStateVectorBodyC::ConvertPointArrayToVec(const Array1dC<Point2dC> &pts) const
{
	UIntT size = pts.Size();
	VectorC res(size*2);
	for(UIntT i = 0; i < size; i ++)
	{
		res[i] = pts[i].Col();
		res[i + size] = pts[i].Row();
	}
	return res;
}

Array1dC<Point2dC> BaseStateVectorBodyC::ProjectAffine(const Array1dC<Point2dC> &pts, const Affine2dC &aff) const
{
	//Project Affine Points
	Array1dC<Point2dC> out(pts.Size());
	for(Array1dIter2C<Point2dC, Point2dC> it(pts,out); it; it++)
	{
		Vector2dC vec(it.Data1().Row(), it.Data1().Col());
		vec = aff * vec;
		Point2dC pt(vec.Row(),vec.Col());
		it.Data2() = pt.Copy();
	}	
	return out;
}


}
