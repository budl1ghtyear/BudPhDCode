#ifndef ICASTATEVECTOR_HH
#define ICASTATEVECTOR_HH

#include "BaseStateVector.hh"
#include "ICAProjection.hh"
#include "SimilarityProjection.hh"
#include "UtilityFunctions.hh"
#include "Ravl/OS/Directory.hh"
namespace RavlN{

class ICAStateVectorBodyC : public BaseStateVectorBodyC
{
	public:
	ICAStateVectorBodyC(){}
	
	ICAStateVectorBodyC(DirectoryC const &dir);
	ICAStateVectorBodyC(DirectoryC const &dir, VectorC const &pts);
	ICAStateVectorBodyC(DirectoryC const &dir, SArray1dC<Point2dC> const &pts);
	ICAStateVectorBodyC(DirectoryC const &dir, Array1dC<Point2dC> const &pts);
	ICAStateVectorBodyC(ICAProjectionC const &p, SimilarityProjectionC const &s, VectorC const &vec):ICA(p),sim(s)
	{this->state = vec.Copy();}

	Array1dC<Point2dC> GetObservationVector();
	
	//Accessor Methods:
	ICAProjectionC GetICA(void) const {return ICA;}
	SimilarityProjectionC GetSimilarityProjection(void) const {return sim;}
	VectorC GetNoiseVector(UIntT const &nwin=1)
	{
		//Affine Noise Portion
		VectorC vec((RealT)nwin,(RealT)nwin,(RealT)1.0,(RealT)1.0);
		return vec.Append(ICA.GetOuterBound(nwin));
	}
	BasisSplineC GetBSpline(void)
	{
		SArray1dC<Point2dC> pts = ArrayToSArray_B(this->GetObservationVector()); //Obtain the B-Spline Control Points  
		BasisSplineC bspl(pts);
		return bspl;
	}
	virtual RCBodyVC &Copy() const
    { return *new ICAStateVectorBodyC(ICA,sim,state.Copy()); }
	protected:
	ICAProjectionC ICA;
	SimilarityProjectionC sim;
	//StringC directory;
};
  
  //! userlevel=Normal
  //: Handle for ICAStateVectorBodyC
  //!cwiz:author
  
  class ICAStateVectorC
    : public BaseStateVectorC
  {
  public:
    ICAStateVectorC() 
      : BaseStateVectorC(*new ICAStateVectorBodyC())
    {}
    //!cwiz:author
    
    ICAStateVectorC(DirectoryC const& dir) 
      : BaseStateVectorC(*new ICAStateVectorBodyC(dir))
    {}
    //!cwiz:author
    
    ICAStateVectorC(DirectoryC const& dir,VectorC const& pts) 
      : BaseStateVectorC(*new ICAStateVectorBodyC(dir,pts))
    {}
    //!cwiz:author
    
    ICAStateVectorC(DirectoryC const& dir,SArray1dC<Point2dC> const& pts) 
      : BaseStateVectorC(*new ICAStateVectorBodyC(dir,pts))
    {}
    //!cwiz:author
    
    ICAStateVectorC(DirectoryC const& dir,Array1dC<Point2dC> const& pts) 
      : BaseStateVectorC(*new ICAStateVectorBodyC(dir,pts))
    {}
    //!cwiz:author
    
    Array1dC<Point2dC> GetObservationVector() 
    { return Body().GetObservationVector(); }
    //Accessor Methods:
    //!cwiz:author
	BasisSplineC GetBSpline()
	{ return Body().GetBSpline();}    
    ICAProjectionC GetICA(void) const
    { return Body().GetICA(); }
    //!cwiz:author
    
    SimilarityProjectionC GetSimilarityProjection(void) const
    { return Body().GetSimilarityProjection(); }
    //!cwiz:author
    
    VectorC GetNoiseVector(UIntT const& nwin = 1) 
    { return Body().GetNoiseVector(nwin); }
    //!cwiz:author
    ICAStateVectorC Copy() const { 
      if(!IsValid()) return ICAStateVectorC();
      return ICAStateVectorC(static_cast<ICAStateVectorBodyC &>(Body().Copy())); 
    }
   
  protected:
    ICAStateVectorC(ICAStateVectorBodyC &bod)
     : BaseStateVectorC(bod)
    {}
    //: Body constructor. 
    
    ICAStateVectorBodyC& Body()
    { return static_cast<ICAStateVectorBodyC &>(BaseStateVectorC::Body()); }
    //: Body Access. 
    
    const ICAStateVectorBodyC& Body() const
    { return static_cast<const ICAStateVectorBodyC &>(BaseStateVectorC::Body()); }
    //: Body Access. 
    
  };
}

#endif
