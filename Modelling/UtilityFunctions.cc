#include "UtilityFunctions.hh"
VectorC PointsToVector(SArray1dC<Point2dC> const &pts)
//Convert points SArray to vector using the new AAMShapeModel notation
// Vect X = R_1,C_1,R_2,C_2....R_n,C_n
{
	VectorC vect(pts.Size()*2);
	for(UIntT i = 0; i < pts.Size(); i++)
	{
		vect[i*2] = pts[i].Row();
		vect[(i*2)+1]=pts[i].Col();
	}
	return vect;
}

VectorC PointsToVector(Array1dC<Point2dC> const &pts)
{
	return PointsToVector(ArrayToSArray_B(pts).Copy());
}

SArray1dC<Point2dC> VectorToPoints(VectorC const &vect)
{
	SArray1dC<Point2dC> pts(vect.Size()/2);
	for(UIntT i = 0; i < pts.Size(); i++)
	{
		Point2dC pt(vect[i*2],vect[(i*2)+1]);
		pts[i] = pt.Copy();
 	}
 	return pts;
}

Point2dC ComputeCentroid(SArray1dC<Point2dC> const &pts)
{
	Point2dC sum(0.0,0.0);
	for(SArray1dIterC<Point2dC> it(pts); it; it++)
		sum += (*it);
	RealT n = (RealT)pts.Size();
	return Point2dC((sum.Row()/n),(sum.Col()/n));
}

