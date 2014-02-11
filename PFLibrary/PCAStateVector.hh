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
	
	PCAStateVectorBodyC(DirectoryC const &dir, UIntT const &numcp = 10); //remember that over here, numcp = number of eigen components to retain in the PCA space
	PCAStateVectorBodyC(DirectoryC const &dir, VectorC const &pts, UIntT const &numcp = 10);
	PCAStateVectorBodyC(DirectoryC const &dir, SArray1dC<Point2dC> const &pts, UIntT const &numcp = 10);//instantiate with the control points of a given B-Spline
	PCAStateVectorBodyC(DirectoryC const &dir, Array1dC<Point2dC> const &pts, UIntT const &numcp = 10);
	PCAStateVectorBodyC(PCAProjectionC const &p, SimilarityProjectionC const &s, VectorC const &vec):pca(p),sim(s)
	{this->state = vec.Copy();}
	Array1dC<Point2dC> GetObservationVector();
	BasisSplineC GetBSpline(void)
	{
		SArray1dC<Point2dC> pts = ArrayToSArray_B(this->GetObservationVector()); //Obtain the B-Spline Control Points  
		BasisSplineC bspl(pts);
		return bspl;
	}
	//Accessor Methods:
	PCAProjectionC GetPCA(void) const {return pca;}
	SimilarityProjectionC GetSimilarityProjection(void) const {return sim;}
	VectorC GetNoiseVector(UIntT const &nwin=1)
	{
		//Affine Noise Portion
		VectorC vec((RealT)nwin,(RealT)nwin,(RealT)1.0,(RealT)1.0);
		return vec.Append(pca.GetOuterBound(nwin));
	}
	   
    virtual RCBodyVC &Copy() const
    { 
		//cout<<"PCASV Copy() Method Called"<<endl;
		return *new PCAStateVectorBodyC(pca.Copy(),sim.Copy(),state.Copy()); 
	}

	protected:
	PCAProjectionC pca;
	SimilarityProjectionC sim;
	//StringC directory;
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
    
 
    Array1dC<Point2dC> GetObservationVector() 
    { return Body().GetObservationVector(); }
	
	BasisSplineC GetBSpline()
	{ return Body().GetBSpline();}
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
    PCAStateVectorC Copy() const { 
      if(!IsValid()) return PCAStateVectorC();
      return PCAStateVectorC(static_cast<PCAStateVectorBodyC &>(Body().Copy())); 
    }
  
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
