#ifndef BASESTATEVECTOR_HH
#define BASESTATEVECTOR_HH
//! author="Bud Goswami"
//! date="22/4/2009"

//Description: This class represents the base state vector
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/Index2d.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/Vector.hh"
namespace RavlN{
  
class BaseStateVectorC
{
public:
	BaseStateVectorC(){} 
	~BaseStateVectorC(){}
	BaseStateVectorC(const VectorC &vec)
	{
		state = vec.Copy();
	}
	//Copy Constructor
	BaseStateVectorC (const BaseStateVectorC &oth)
	{
		state = oth.state.Copy();
	}
	virtual BaseStateVectorC& operator=(const BaseStateVectorC &oth)
	{
		state = oth.state;
		return *this;
	}
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
		return obs;
	}
	//This Virtual Function Ensures that the different state vectors have individual obs vectors
	virtual VectorC GetStateVector() {return state;} 
	//This Virtual Function Ensures that different State Vector Classes can implement their own State vectors

protected:
	VectorC state;	
};
  

}

#endif

