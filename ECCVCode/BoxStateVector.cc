#include "BoxStateVector.hh"
namespace RavlN
{
	Array1dC<Point2dC> BoxStateVectorBodyC::GetObservationVector()
	{
		Index2dC center(state[0],state[1]);
		IndexRange2dC range(center, numrows,numrows);
		Array1dC<Point2dC> res(4);
		res[0] = Point2dC(range.TopLeft().Row(),range.TopLeft().Col()).Copy();
		res[1] = Point2dC(range.TopRight().Row(),range.TopRight().Col()).Copy();
		res[2] = Point2dC(range.BottomLeft().Row(),range.BottomLeft().Col()).Copy();
		res[3] = Point2dC(range.BottomRight().Row(),range.BottomRight().Col()).Copy();
		return res;
	}
}
