#ifndef EIGENSTATEVECTOR_HH
#define EIGENSTATEVECTOR_HH
//! author="Bud Goswami"
//! date="20/8/2008"
//This class represents an eigen-state vector
//An external, head affine transform is used to geometrically normalise the lip co-ordinates.
//Subsequent projection onto the provided eigen space parameters leads to the determination of the state space
#include "BaseStateVector.hh"

namespace RavlN {
  
class EigenStateVectorBodyC : public BaseStateVectorBodyC
{
public:
	EigenStateVectorBodyC() {} //Default Constructor
	EigenStateVectorBodyC(const Array1dC<Point2dC> &pts, const Affine2dC &headaff,MatrixC &mn, MatrixC &ev, const UIntT &dim);
	EigenStateVectorBodyC(const Array1dC<Point2dC> &pts, const Affine2dC &headaff,MatrixC &mn, MatrixC &ev);
	EigenStateVectorBodyC(const Affine2dC &headaff,const MatrixC &mn,const MatrixC &ev,const MatrixC &evl, const UIntT &nm, const VectorC &v): BaseStateVectorBodyC(v.Copy())
	{
		head_aff = headaff;
		mean = mn;
		evect = ev;
		evals = evl.Copy();
		numeig = nm;
	}	
	virtual RCBodyVC &Copy() const
   { return *new EigenStateVectorBodyC(head_aff,mean,evect,evals.Copy(),numeig,state.Copy()); }

	//The standard state vector methods	
	Array1dC<Point2dC> GetPreAffineObservationVector();
	Array1dC<Point2dC> GetObservationVector();
	//VectorC GetStateVector() {return state;}
	//Accessor Methods
	MatrixC GetEigenValues(void) const {return evals;}
	Affine2dC GetHeadAffineProjection(void) const {return head_aff;}
	void SetHeadAffine(const Affine2dC &aff) {head_aff = aff;}
	UIntT GetNumEigDimensions(void) const {return numeig;}
	void SetNumEigDimensions(const UIntT &n) {numeig = n;} 
	MatrixC GetMeanMatrix(void) const {return mean;}
	MatrixC GetEigenVectors(void) const {return evect;}

protected:
	
	MatrixC mean;
	MatrixC evect;
	MatrixC evals;
	Affine2dC head_aff;
	UIntT numeig;
	MatrixC ComputeEigenValues(const VectorC &pt_vec, const MatrixC &mean, const MatrixC &evect);

};
  
  //! userlevel=Normal
  //: Handle for EigenStateVectorBodyC
  //!cwiz:author
  
  class EigenStateVectorC
    : public BaseStateVectorC
  {
  public:
    EigenStateVectorC() 
      : BaseStateVectorC(*new EigenStateVectorBodyC())
    {}
    //Default Constructor
    //!cwiz:author
    
    EigenStateVectorC(const Array1dC<Point2dC> & pts,const Affine2dC & headaff,MatrixC & mn,MatrixC & ev,const UIntT & dim) 
      : BaseStateVectorC(*new EigenStateVectorBodyC(pts,headaff,mn,ev,dim))
    {}
    //!cwiz:author
    
    EigenStateVectorC(const Array1dC<Point2dC> & pts,const Affine2dC & headaff,MatrixC & mn,MatrixC & ev) 
      : BaseStateVectorC(*new EigenStateVectorBodyC(pts,headaff,mn,ev))
    {}
    //!cwiz:author
    
    EigenStateVectorC(const Affine2dC & headaff,const MatrixC & mn,const MatrixC & ev,const MatrixC & evl,const UIntT & nm,const VectorC & v) 
      : BaseStateVectorC(*new EigenStateVectorBodyC(headaff,mn,ev,evl,nm,v))
    {}
    //!cwiz:author
	 EigenStateVectorC Copy() const 
	 { 
      if(!IsValid()) return EigenStateVectorC();
      return EigenStateVectorC(static_cast<EigenStateVectorBodyC &>(Body().Copy())); 
    }
    Array1dC<Point2dC> GetPreAffineObservationVector() 
    { return Body().GetPreAffineObservationVector(); }
    //!cwiz:author
    
    Array1dC<Point2dC> GetObservationVector() 
    { return Body().GetObservationVector(); }
    //VectorC GetStateVector() {return state;}
    //Accessor Methods
    //!cwiz:author
    
    MatrixC GetEigenValues(void) const
    { return Body().GetEigenValues(); }
    //!cwiz:author
    
    Affine2dC GetHeadAffineProjection(void) const
    { return Body().GetHeadAffineProjection(); }
    //!cwiz:author
    
    void SetHeadAffine(const Affine2dC & aff) 
    { Body().SetHeadAffine(aff); }
    //!cwiz:author
    
    UIntT GetNumEigDimensions(void) const
    { return Body().GetNumEigDimensions(); }
    //!cwiz:author
    
    void SetNumEigDimensions(const UIntT & n) 
    { Body().SetNumEigDimensions(n); }
    //!cwiz:author
    
    MatrixC GetMeanMatrix(void) const
    { return Body().GetMeanMatrix(); }
    //!cwiz:author
    
    MatrixC GetEigenVectors(void) const
    { return Body().GetEigenVectors(); }
    //Overloaded operators
    //!cwiz:author
    
    //Copy Constructor
    //!cwiz:author
    
    EigenStateVectorBodyC & operator=(const EigenStateVectorBodyC & oth) 
    { return Body().operator=(oth); }
    //Copy Constructor
    //!cwiz:author
    
  protected:
    EigenStateVectorC(EigenStateVectorBodyC &bod)
     : BaseStateVectorC(bod)
    {}
    //: Body constructor. 
    
    EigenStateVectorBodyC& Body()
    { return static_cast<EigenStateVectorBodyC &>(BaseStateVectorC::Body()); }
    //: Body Access. 
    
    const EigenStateVectorBodyC& Body() const
    { return static_cast<const EigenStateVectorBodyC &>(BaseStateVectorC::Body()); }
    //: Body Access. 
    
  };


}

#endif
