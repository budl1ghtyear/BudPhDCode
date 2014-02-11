#ifndef BOXSTATEVECTOR_HH
#define BOXSTATEVECTOR_HH

#include "BaseStateVector.hh"
#include "Ravl/OS/Directory.hh"
namespace RavlN
{

class BoxStateVectorBodyC : public BaseStateVectorBodyC
{
	public:
	BoxStateVectorBodyC(){;}
	BoxStateVectorBodyC(Index2dC const &ctr, UIntT const &rows=102, UIntT const &cols=122):numrows(rows),numcols(cols)
	{
		VectorC st(2);st.Fill(0.);
		st[0] = ctr.Row();
		st[1] = ctr.Col();
		state = st.Copy();
	}
	BoxStateVectorBodyC(VectorC const &st, UIntT const &rows=102, UIntT const &cols=122):BaseStateVectorBodyC(st.Copy()),numrows(rows),numcols(cols)
	{}
	Array1dC<Point2dC> GetObservationVector();
	//Accessor Methods:
	BasisSplineC GetBSpline()
	{
		BasisSplineC spl(11,2,30);
		return spl;
	}
	VectorC GetNoiseVector(UIntT const &nwin=5)
	{
		//Affine Noise Portion
		VectorC vec((RealT)nwin,(RealT)nwin);
		return vec;
	}
    virtual RCBodyVC &Copy() const
    { 
		//cout<<"PCASV Copy() Method Called"<<endl;
		return *new BoxStateVectorBodyC(state.Copy(),numrows.V(),numcols.V()); 
	}
	void SetNumRows(IndexC const &nr) {numrows = nr;}
	void SetNumCols(IndexC const &nr) {numcols = nr;}
	protected:
	IndexC numrows;
	IndexC numcols;
	//StringC directory;
};
  
  //! userlevel=Normal
  //: Handle for BoxStateVectorBodyC
  //!cwiz:author
  
  class BoxStateVectorC
    : public BaseStateVectorC
  {
  public:
    BoxStateVectorC() 
      : BaseStateVectorC(*new BoxStateVectorBodyC())
    {}
    //!cwiz:author
    
    BoxStateVectorC(Index2dC const& ctr,UIntT const& rows = 102,UIntT const& cols = 122) 
      : BaseStateVectorC(*new BoxStateVectorBodyC(ctr,rows,cols))
    {}
    //!cwiz:author
    
    BoxStateVectorC(VectorC const& st,UIntT const& rows = 102,UIntT const& cols = 122) 
      : BaseStateVectorC(*new BoxStateVectorBodyC(st,rows,cols))
    {}
    //!cwiz:author
    
    Array1dC<Point2dC> GetObservationVector() 
    { return Body().GetObservationVector(); }
    //Accessor Methods:
    //!cwiz:author
    
    BasisSplineC GetBSpline() 
    { return Body().GetBSpline(); }
    //!cwiz:author
    
    VectorC GetNoiseVector(UIntT const& nwin = 5) 
    { return Body().GetNoiseVector(nwin); }
    //!cwiz:author
    
    RCBodyVC & Copy() const
    { return Body().Copy(); }
    //!cwiz:author
    
    void SetNumRows(IndexC const& nr) 
    { Body().SetNumRows(nr); }
    //!cwiz:author
    
    void SetNumCols(IndexC const& nr) 
    { Body().SetNumCols(nr); }
    //!cwiz:author
    
  protected:
    BoxStateVectorC(BoxStateVectorBodyC &bod)
     : BaseStateVectorC(bod)
    {}
    //: Body constructor. 
    
    BoxStateVectorBodyC& Body()
    { return static_cast<BoxStateVectorBodyC &>(BaseStateVectorC::Body()); }
    //: Body Access. 
    
    const BoxStateVectorBodyC& Body() const
    { return static_cast<const BoxStateVectorBodyC &>(BaseStateVectorC::Body()); }
    //: Body Access. 
    
  };


  
  
}

#endif
