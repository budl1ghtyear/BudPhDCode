#include "PointStateVectorC.hh"
namespace RavlN{

Array1dC<Point2dC> PointStateVectorBodyC::GetPointsArray(void)
{
	UIntT size = this->x.Size()/2;
	Array1dC<Point2dC> pts(size);
	for(UIntT i = 0; i < size; i++)
	{
		Point2dC pt(0,0);
		pt.Col() = this->x[i];
		pt.Row() = this->x[i + size];
		pts[i] = pt.Copy();
	}	
	return pts;
}

}

