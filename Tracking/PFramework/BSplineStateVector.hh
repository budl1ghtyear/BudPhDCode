#ifndef BSPLINESTATEVECTOR_H
#define BSPLINESTATEVECTOR_H

#include "BaseStateVector.hh"
#include "BSplineC.hh"
namespace RavlN{

class BSplineStateVectorBodyC : public BaseStateVectorBodyC, BSplineC
{
	public:
	BSplineStateVectorBodyC(){}
	virtual ~BSplineStateVectorBodyC(){}
	
	BSplineStateVectorBodyC(const Array1dC<Point2dC> &control_pts, const Array1dC<RealT> &knot_vec, const UIntT &num_pts, const UIntT &ord);
	Array1dC<Point2dC> GetObservationVector();
	BSplineC GetBSpline() 
	{
		return BSplineC(BSplineC::order, BSplineC::knot_vec, BSplineC::ctrl_vec);
	}
	virtual RCBodyVC &Copy() const
    { return *new BSplineStateVectorBodyC(BaseStateVectorBodyC::ConvertVecToPointArray(this->state.Copy()),this->BSplineC::GetKnotVector(),num_points,this->BSplineC::GetOrder());}
	protected:
	//Array1dC<Point2dC> ctrl_points;	
	//Array1dC<RealT> k_vector;	
	UIntT num_points;
	//UIntT order;
	//BSplineC bSpline;
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
  //: Handle for BaseStateVectorBodyC
  //!cwiz:author
  
  class BSplineStateVectorC
    : public BaseStateVectorC
  {
  public:
    BSplineStateVectorC() : BaseStateVectorC(*new BSplineStateVectorBodyC())
    {}
    //!cwiz:author
	BSplineStateVectorC(const Array1dC<Point2dC> &control_pts, const Array1dC<RealT> &knot_vec, const UIntT &num_pts = 20, const UIntT &ord = 4) : 
	BaseStateVectorC(*new BSplineStateVectorBodyC(control_pts, knot_vec, num_pts, ord)) {}
	   
	BSplineStateVectorC Copy() const 
	{ 
      if(!IsValid()) return BSplineStateVectorC();
      return BSplineStateVectorC(static_cast<BSplineStateVectorBodyC &>(Body().Copy())); 
    }
    
	Array1dC<Point2dC> GetObservationVector() 
    { return Body().GetObservationVector(); }
	BSplineC GetBSpline() 
	{
		return Body().GetBSpline();
	}
  protected:
    BSplineStateVectorC(BSplineStateVectorBodyC &bod)
     : BaseStateVectorC(bod)
    {}
    //: Body constructor. 
    
    BSplineStateVectorBodyC& Body()
    { return static_cast<BSplineStateVectorBodyC &>(BaseStateVectorC::Body()); }
    //: Body Access. 
    
    const BSplineStateVectorBodyC& Body() const
    { return static_cast<const BSplineStateVectorBodyC &>(BaseStateVectorC::Body()); }
    //: Body Access. 

  };
}
#endif
