#ifndef PCASTATEVECTOR_HH
#define PCASTATEVECTOR_HH

#include "BaseStateVector.hh"
#include "PCAProjection.hh"
#include "SimilarityProjection.hh"
#include "UtilityFunctions.hh"
#include "Ravl/OS/Directory.hh"
namespace RavlN{

class PCAStateVectorBodyC : public BaseStateVectorBodyC
{
	public:
	PCAStateVectorBodyC(){}
	
	PCAStateVectorBodyC(DirectoryC const &dir, UIntT const &numcp = 10);
	PCAStateVectorBodyC(DirectoryC const &dir, VectorC const &pts, UIntT const &numcp = 10);
	PCAStateVectorBodyC(DirectoryC const &dir, SArray1dC<Point2dC> const &pts, UIntT const &numcp = 10);
	PCAStateVectorBodyC(DirectoryC const &dir, Array1dC<Point2dC> const &pts, UIntT const &numcp = 10);
	virtual RCBodyVC &Copy() const
    { return *new PCAStateVectorBodyC(dir,state.Copy()); }
	Array1dC<Point2dC> GetObservationVector();
	
	//Accessor Methods:
	PCAProjectionC GetPCA(void) const {return pca;}
	SimilarityProjectionC GetSimilarityProjection(void) const {return sim;}
	VectorC GetNoiseVector(UIntT const &nwin=1)
	{
		VectorC vec((RealT)nwin,(RealT)nwin,(RealT)nwin,(RealT)nwin);
		return vec.Append(pca.GetOuterBound(nwin));
	}
	protected:
	PCAProjectionC pca;
	SimilarityProjectionC sim;
	DirectoryC dir;
};
  
  //! userlevel=Normal
  //: Handle for PCAStateVectorBodyC
  //!cwiz:author
  
  class PCAStateVectorC
    : public BaseStateVectorC
  {
  public:
    PCAStateVectorC() 
      : BaseStateVectorC(*new PCAStateVectorBodyC())
    {}
    //!cwiz:author
    
    PCAStateVectorC(DirectoryC const& dir,UIntT const& numcp = 10) 
      : BaseStateVectorC(*new PCAStateVectorBodyC(dir,numcp))
    {}
    //!cwiz:author
    
    PCAStateVectorC(DirectoryC const& dir,VectorC const& pts,UIntT const& numcp = 10) 
      : BaseStateVectorC(*new PCAStateVectorBodyC(dir,pts,numcp))
    {}
    //!cwiz:author
    
    PCAStateVectorC(DirectoryC const& dir,SArray1dC<Point2dC> const& pts,UIntT const& numcp = 10) 
      : BaseStateVectorC(*new PCAStateVectorBodyC(dir,pts,numcp))
    {}
    //!cwiz:author
    
    PCAStateVectorC(DirectoryC const& dir,Array1dC<Point2dC> const& pts,UIntT const& numcp = 10) 
      : BaseStateVectorC(*new PCAStateVectorBodyC(dir,pts,numcp))
    {}
    //!cwiz:author
    
    RCBodyVC & Copy() const
    { return Body().Copy(); }
    //!cwiz:author
    
    Array1dC<Point2dC> GetObservationVector() 
    { return Body().GetObservationVector(); }
    //Accessor Methods:
    //!cwiz:author
    
    PCAProjectionC GetPCA(void) const
    { return Body().GetPCA(); }
    //!cwiz:author
    
    SimilarityProjectionC GetSimilarityProjection(void) const
    { return Body().GetSimilarityProjection(); }
    //!cwiz:author
    
    VectorC GetNoiseVector(UIntT const& nwin = 1) 
    { return Body().GetNoiseVector(nwin); }
    //!cwiz:author
    
  protected:
    PCAStateVectorC(PCAStateVectorBodyC &bod)
     : BaseStateVectorC(bod)
    {}
    //: Body constructor. 
    
    PCAStateVectorBodyC& Body()
    { return static_cast<PCAStateVectorBodyC &>(BaseStateVectorC::Body()); }
    //: Body Access. 
    
    const PCAStateVectorBodyC& Body() const
    { return static_cast<const PCAStateVectorBodyC &>(BaseStateVectorC::Body()); }
    //: Body Access. 
    
  };
}

#endif
