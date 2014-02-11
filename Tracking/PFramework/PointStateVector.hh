#ifndef POINTSTATEVECTOR_HH
#define POINTSTATEVECTOR_HH
//! author="Bud Goswami"
//! date="20/8/2008"

#include "BaseStateVector.hh"

namespace RavlN {
  
class PointStateVectorBodyC : public BaseStateVectorBodyC
{
public:
	PointStateVectorBodyC(){}
	
	PointStateVectorBodyC(Array1dC<Point2dC> &ar):BaseStateVectorBodyC(ConvertPointToVec(ar)){}
	
	PointStateVectorBodyC(Array1dC<Index2dC> &ar):BaseStateVectorBodyC(ConvertPointToVec(ar)){}
	PointStateVectorBodyC(const VectorC &v):BaseStateVectorBodyC(v)
	{}
	
	virtual RCBodyVC &Copy() const
   { return *new PointStateVectorBodyC(state.Copy());}
   
	virtual Array1dC<Point2dC> GetObservationVector()
	{	
		//Assumes that the VectorC is the same as the Point configuration
		UIntT size = state.Size()/2;
		Array1dC<Point2dC> obs(size);
		for(UIntT i = 0; i < size; i++)
		{
			Point2dC pt(0.0,0.0);
			pt.Col() = state[i];
			pt.Row() = state[i + size];
			obs[i] = pt.Copy();
		}
		//Point2dC pt(10,10);
		//Array1dC<Point2dC> obs(1); obs[0] = pt;
		//cout<<"Calling BaseStateVectorC::GetObservationVector()"<<endl;
		return obs;
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
  
  //! userlevel=Normal
  //: Handle for PointStateVectorBodyC
  //!cwiz:author
  
  class PointStateVectorC
    : public BaseStateVectorC
  {
  public:
    PointStateVectorC() 
      : BaseStateVectorC(*new PointStateVectorBodyC())
    {}
    //!cwiz:author
    
    PointStateVectorC(Array1dC<Point2dC> & ar) 
      : BaseStateVectorC(*new PointStateVectorBodyC(ar))
    {}
    //!cwiz:author
    
    PointStateVectorC(Array1dC<Index2dC> & ar) 
      : BaseStateVectorC(*new PointStateVectorBodyC(ar))
    {}
    //!cwiz:author
    
    PointStateVectorC(const VectorC & v) 
      : BaseStateVectorC(*new PointStateVectorBodyC(v))
    {}
    //!cwiz:author
    
	 PointStateVectorC Copy() const 
	 { 
      if(!IsValid()) return PointStateVectorC();
      return PointStateVectorC(static_cast<PointStateVectorBodyC &>(Body().Copy())); 
    }
    
    Array1dC<Point2dC> GetObservationVector() 
    { return Body().GetObservationVector(); }
    //!cwiz:author
    
    Array1dC<Point2dC> GetPointsArray(void) 
    { return Body().GetPointsArray(); }
    //!cwiz:author
    
  protected:
    PointStateVectorC(PointStateVectorBodyC &bod)
     : BaseStateVectorC(bod)
    {}
    //: Body constructor. 
    
    PointStateVectorBodyC& Body()
    { return static_cast<PointStateVectorBodyC &>(BaseStateVectorC::Body()); }
    //: Body Access. 
    
    const PointStateVectorBodyC& Body() const
    { return static_cast<const PointStateVectorBodyC &>(BaseStateVectorC::Body()); }
    //: Body Access. 
    
  };


}

#endif
