#ifndef EASTATEVECTOR_HH
#define EASTATEVECTOR_HH

#include "BSplineStateVector.hh"
#include 	"Ravl/PatternRec/FuncMeanProjection.hh"
namespace RavlN{

class EAStateVectorBodyC : public BaseStateVectorBodyC, BSplineC//BSplineStateVectorBodyC
{
	public:
	EAStateVectorBodyC(){}
	virtual ~EAStateVectorBodyC() {}
	
	EAStateVectorBodyC(const Array1dC<Point2dC> &mean, const MatrixC &evect, const Affine2dC &aff, const Array1dC<Point2dC> &ctrl_pts, const Array1dC<RealT> &knot_vec, const UIntT &ord = 4, const UIntT &n_pts = 20);
	
	Array1dC<Point2dC> GetObservationVector();
	virtual RCBodyVC &Copy() const
    { return *new EAStateVectorBodyC(BaseStateVectorBodyC::ConvertVecToPointArray(this->pca.Mean()).Copy(), (MatrixC)pca.Projection().Copy(), this->GetAffineTransform(), this->GetControlPoints().Copy(), this->BSplineC::GetKnotVector(),this->BSplineC::GetOrder(), num_points);}
	
	Array1dC<Point2dC> GetControlPoints() const;
	Affine2dC GetAffineTransform() const;
	BSplineC GetBSpline()
	{
		return BSplineC(BSplineC::order, BSplineC::knot_vec, BSplineC::ctrl_vec);
	}
	protected:
	UIntT n_pts;
	FuncMeanProjectionC pca;
	UIntT num_points;
	//static VectorC ConvertPointToVec(const Array1dC<Point2dC> &pts, const Array1dC<Point2dC> &mean, const MatrixC &evect, const Affine2dC &af)
	VectorC ConvertPointToVec(const Array1dC<Point2dC> &pts, const Array1dC<Point2dC> &mean, const MatrixC &evect, const Affine2dC &af)
	{
		//Project the points using the affine transform
		Array1dC<Point2dC> c_pts = pts.Copy();
		for(Array1dIterC<Point2dC> it(c_pts); it; it++)
		{
			Vector2dC vec((*it).Row(),(*it).Col());
			vec = af * vec;
			(*it) = (Point2dC(vec.Row(),vec.Col())).Copy();
		}
		Point2dC translation;Vector2dC scaling;RealT skew;RealT rotation;
		af.Decompose(translation,scaling, skew, rotation);
		VectorC af_components (translation.Row(),translation.Col(),scaling.Row(),scaling.Col(),skew,rotation);
		//Now convert to a vector and obtain the eigenvalues
		//Convert the array of points to a vector
		//The following is the most beautiful use of RAVL I have encountered thus far! Its how prototyping SHOULD be allowed!
		//cout<<"INSIDE EASTATEVECTOR - Affine transformed cpts are: \n" << c_pts<<endl;
		//cout<<"INSIDE EASTATEVECTOR - Affine transform is: \n" << af<<endl;
		VectorC cpt_vector = BaseStateVectorBodyC::ConvertPointArrayToVec(pts);
		VectorC mean_vect = BaseStateVectorBodyC::ConvertPointArrayToVec(mean);
		FuncMeanProjectionC pca(mean_vect, evect);
		//VectorC evals = pca.Apply(cpt_vector);
		VectorC evals = evect.Inverse() * (cpt_vector - mean_vect);
		
		return evals.Append(af_components);
	}
};

  class EAStateVectorC
    : public BaseStateVectorC
  {
  public:
    EAStateVectorC() : BaseStateVectorC(*new EAStateVectorBodyC())
    {}
    //!cwiz:author
	EAStateVectorC(const Array1dC<Point2dC> &mean, const MatrixC &evect, const Affine2dC &aff, const Array1dC<Point2dC> &ctrl_pts, const Array1dC<RealT> &knot_vec, const UIntT &ord = 4, const UIntT &n_pts = 20) :
	BaseStateVectorC(*new EAStateVectorBodyC(mean,evect,aff,ctrl_pts,knot_vec,ord,n_pts)) {}
	   
	EAStateVectorC Copy() const 
	{ 
      if(!IsValid()) return EAStateVectorC();
      return EAStateVectorC(static_cast<EAStateVectorBodyC &>(Body().Copy())); 
    }
    
    Array1dC<Point2dC> GetControlPoints() const
    { return Body().GetControlPoints(); }
    //!cwiz:author
    
    Affine2dC GetAffineTransform() const
    { return Body().GetAffineTransform(); }
    //!cwiz:author
    
	Array1dC<Point2dC> GetObservationVector() 
    { return Body().GetObservationVector(); }
	
	Array1dC<Point2dC> GetControlPoints()
	{return Body().GetControlPoints();}
	
	Affine2dC GetAffineTransform()
	{return Body().GetAffineTransform();}
	BSplineC GetBSpline() 
	{
		return Body().GetBSpline();
	}
	
  protected:
    EAStateVectorC(EAStateVectorBodyC &bod)
     : BaseStateVectorC(bod)
    {}
    //: Body constructor. 
    
    EAStateVectorBodyC& Body()
    { return static_cast<EAStateVectorBodyC &>(BaseStateVectorC::Body()); }
    //: Body Access. 
    
    const EAStateVectorBodyC& Body() const
    { return static_cast<const EAStateVectorBodyC &>(BaseStateVectorC::Body()); }
    //: Body Access. 

  };
	
}

#endif
