#ifndef ABSTRACTINH_HH
#define ABSTRACTINH_HH
//! author="Bud Goswami"
//! date="29/5/2009"

//Description: This class represents the abstract base class
#include "Ravl/RefCounter.hh"

namespace RavlN{
  
class AbstractInhBodyC : public RCBodyVC

{
public:
	AbstractInhBodyC(){} 
	~AbstractInhBodyC(){}

	virtual RCBodyVC &Copy() const = 0;
	virtual IntT AbstractFunction() = 0;
};
  
  //! userlevel=Normal
  //: Handle for AbstractInhBodyC
  //!cwiz:author
  
  class AbstractInhC
    : public RCHandleC<AbstractInhBodyC>
  {
  public:
    AbstractInhC() 
      : RCHandleC<AbstractInhBodyC>(*new AbstractInhBodyC())
    {}
    //!cwiz:author
    //The above method is WRONG. A pure virtual function cannot have an object or a pointer instantiation statement. In the example, the following is suggested:
    //AbstractInhC(){} // Invalid default constructor
    
    RCBodyVC & Copy() const
    { return Body().Copy(); }
    //!cwiz:author
    //Also causes issues. NEEDS to be as below:
    //BaseStateVectorC Copy() const { 
    //  if(!IsValid()) return BaseStateVectorC();
    //  return BaseStateVectorC(static_cast<BaseStateVectorBodyC &>(Body().Copy())); 
    //  }
    
    IntT AbstractFunction() 
    { return Body().AbstractFunction(); }
    //!cwiz:author
    
  protected:
    AbstractInhC(AbstractInhBodyC &bod)
     : RCHandleC<AbstractInhBodyC>(bod)
    {}
    //: Body constructor. 
    
    AbstractInhBodyC& Body()
    { return static_cast<AbstractInhBodyC &>(RCHandleC<AbstractInhBodyC>::Body()); }
    //: Body Access. 
    
    const AbstractInhBodyC& Body() const
    { return static_cast<const AbstractInhBodyC &>(RCHandleC<AbstractInhBodyC>::Body()); }
    //: Body Access. 
    
  };

}
#endif

