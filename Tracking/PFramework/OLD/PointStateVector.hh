#ifndef POINTSTATEVECTOR_HH
#define POINTSTATEVECTOR_HH
//! author="Bud Goswami"
//! date="20/8/2008"

#include "BaseStateVector.hh"
#include "Ravl/Matrix3d.hh"
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/Index2d.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/Vector.hh"

namespace RavlN {
  
class PointStateVectorC : public BaseStateVectorC
{
public:
	PointStateVectorC(){}
	
	PointStateVectorC(Array1dC<Point2dC> &ar):BaseStateVectorC(ConvertPointToVec(ar)){}
	
	PointStateVectorC(Array1dC<Index2dC> &ar):BaseStateVectorC(ConvertPointToVec(ar)){}
	
	/*PointStateVectorC& operator=(const PointStateVectorC &oth)
	{
		state = oth.state;
		return *this;
	}*/
	
	PointStateVectorC(const PointStateVectorC &pt)
	{
		state = pt.state.Copy();
	}
	
	Array1dC<Point2dC> GetPointsArray(void);
protected:
	static VectorC ConvertPointToVec(Array1dC<Point2dC> ar)
	{
		//This method converts to points to a vector {Xo,...,Xm,Yo,...Ym}
		//m = total number of points
		//X = horizontal dimension, Cols(RAVL Indexing)
		//Y = vertical dimension, Rows(RAVL Indexing)
		VectorC out((ar.Size()*2));
		out.Fill(0.0);
		IndexC i = 0;
		for(Array1dIterC<Point2dC> it(ar); it; it++)
		{
			out[i] = (RealT)((*it).Col());
			out[i+ar.Size()]= (RealT)((*it).Row());
			i++;
		}
		return out;
	}
	static VectorC ConvertPointToVec(Array1dC<Index2dC> ar)
	{
		//This method converts to points to a vector {Xo,...,Xm,Yo,...Ym}
		//m = total number of points
		//X = horizontal dimension, Cols(RAVL Indexing)
		//Y = vertical dimension, Rows(RAVL Indexing)
		VectorC out((ar.Size()*2));
		out.Fill(0.0);
		IndexC i = 0;
		for(Array1dIterC<Index2dC> it(ar); it; it++)
		{
			out[i] = (RealT)((*it).Col());
			out[i+ar.Size()]= (RealT)((*it).Row());
			i++;
		}
		return out;
	}
};


}

#endif
