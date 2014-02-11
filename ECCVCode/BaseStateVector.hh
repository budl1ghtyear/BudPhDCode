#ifndef BASESTATEVECTOR_HH
#define BASESTATEVECTOR_HH
//! author="Bud Goswami"
//! date="22/4/2009"

//Description: This class represents the base state vector
#include "Ravl/Affine2d.hh"
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/Array1dIter2.hh"
#include "Ravl/Matrix.hh"
#include "Ravl/Index2d.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/RefCounter.hh"
#include "Ravl/Vector.hh"
//#include "BSplineC.hh"
#include "BasisSpline.hh"
#include "UtilityFunctions.hh"
namespace RavlN{
  
class BaseStateVectorBodyC : public RCBodyVC

{
public:
	BaseStateVectorBodyC(){} 
	virtual ~BaseStateVectorBodyC(){}
	BaseStateVectorBodyC(const VectorC &vec)
	{
		state = vec.Copy();
	}


	virtual Array1dC<Point2dC> GetObservationVector() = 0;
	//virtual BSplineC GetBSpline() = 0;
	//This Virtual Function Ensures that the different state vectors have individual obs vectors
	VectorC GetStateVector() {return state;} 
	virtual void SetStateVector(const VectorC &vec) 
	{
		if((vec.Size() == state.Size())&&(vec.IsEmpty()!=true)) 
		state = vec;
	}
	virtual VectorC GetNoiseVector(UIntT const &nwin=1) = 0;
	virtual BasisSplineC GetBSpline(void) = 0;
	//virtual operator=
	//This Virtual Function Ensures that different State Vector Classes can implement their own State vectors
	//~ virtual BaseStateVectorBodyC& operator=(const BaseStateVectorBodyC &oth) = 0;
	//~ {	
		//~ if(this != &oth)
		//~ {
		//~ state = oth.state.Copy();
		//~ }
		//~ return *this;
	//~ }
	virtual RCBodyVC &Copy() const = 0;

protected:
	VectorC state;	
	Array1dC<Point2dC> ConvertVecToPointArray(const VectorC &v) const;
	VectorC ConvertPointArrayToVec(const Array1dC<Point2dC> &pts) const;
	Array1dC<Point2dC> ProjectAffine(const Array1dC<Point2dC> &pts, const Affine2dC &aff) const;
};
  
  //! userlevel=Normal
  //: Handle for BaseStateVectorBodyC
  //!cwiz:author
  
  class BaseStateVectorC
    : public RCHandleC<BaseStateVectorBodyC>
  {
  public:
    BaseStateVectorC() 
    {}
    //!cwiz:author
    
    BaseStateVectorC(const VectorC & vec) 
    {}
    //!cwiz:author
	
    Array1dC<Point2dC> GetObservationVector() 
    { return Body().GetObservationVector(); }
	BasisSplineC GetBSpline() 
    { return Body().GetBSpline(); }
    //This Virtual Function Ensures that the different state vectors have individual obs vectors
    //!cwiz:author
    //~ BSplineC GetBSpline()
	//~ {return Body().GetBSpline();}

    VectorC GetStateVector() 
    { return Body().GetStateVector(); }
    //!cwiz:author
    
    virtual void SetStateVector(const VectorC & vec) 
    { Body().SetStateVector(vec); }
    //This Virtual Function Ensures that different State Vector Classes can implement their own State vectors
    //!cwiz:author
    
    VectorC GetNoiseVector(UIntT const& nwin = 1) 
    { return Body().GetNoiseVector(nwin); }
    //virtual operator=
    //This Virtual Function Ensures that different State Vector Classes can implement their own State vectors
    //virtual BaseStateVectorBodyC& operator=(const BaseStateVectorBodyC &oth) = 0;
    //!cwiz:author
    
    BaseStateVectorBodyC & operator=(const BaseStateVectorBodyC & oth) 
    { return Body().operator=(oth); }
    //!cwiz:author
    BaseStateVectorC Copy() const { 
      if(!IsValid()) return BaseStateVectorC();
      return BaseStateVectorC(static_cast<BaseStateVectorBodyC &>(Body().Copy())); 
    }
   
  protected:
    BaseStateVectorC(BaseStateVectorBodyC &bod)
     : RCHandleC<BaseStateVectorBodyC>(bod)
    {}
    //: Body constructor. 
    
    BaseStateVectorBodyC& Body()
    { return static_cast<BaseStateVectorBodyC &>(RCHandleC<BaseStateVectorBodyC>::Body()); }
    //: Body Access. 
    
    const BaseStateVectorBodyC& Body() const
    { return static_cast<const BaseStateVectorBodyC &>(RCHandleC<BaseStateVectorBodyC>::Body()); }
    //: Body Access. 
    
  };
  

}

#endif

