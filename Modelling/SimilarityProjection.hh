#ifndef SimilarityPROJECTION_HH
#define SimilarityPROJECTION_HH

//class to perform affine projection and inverse projection

#include "UtilityFunctions.hh"
#include "Ravl/Affine2d.hh"
#include "Ravl/Vector2d.hh"
#include "Ravl/Tuple2.hh"
using namespace RavlN;

class SimilarityProjectionC
{
	public:
	SimilarityProjectionC(){}
	~SimilarityProjectionC(){}
	
	SimilarityProjectionC(StringC const &dir);
	SimilarityProjectionC(SArray1dC<Point2dC> const &affinemean);
	SimilarityProjectionC(StringC const &str, SArray1dC<Point2dC> const &am, Affine2dC const &aff):directory(str),a_mean(am),sim_proj(aff) {;}
	SimilarityProjectionC(SimilarityProjectionC const &proj)
	{
		directory = proj.GetDirectory();
		a_mean = proj.GetAffineMean();
		sim_proj = proj.GetAffine();
		//cout<<"SimilarityProjectionC Copy Called"<<endl;
	}
	SimilarityProjectionC Copy() const
	{
		return SimilarityProjectionC(directory,a_mean, sim_proj);
	}
	Tuple2C<SArray1dC<Point2dC>,VectorC> ProjectToAffineSpace(SArray1dC<Point2dC> const &pts); //perform a similarity projection onto the mean and returns the result of the residual of the projection
	SArray1dC<Point2dC> ProjectFromAffineSpace(VectorC const &vec, VectorC const &aff);
	SArray1dC<Point2dC> ProjectFromAffineSpace(SArray1dC<Point2dC> const &vec, VectorC const &aff)
	{return this->ProjectFromAffineSpace(PointsToVector(vec),aff);}
	SArray1dC<Point2dC> GetMean(void) const {return a_mean;}
	//Accessor methods:
	StringC GetDirectory(void) const {return directory;}
	SArray1dC<Point2dC> GetAffineMean(void) const {return a_mean;}
	Affine2dC GetAffine(void) const {return sim_proj;}
		
	protected:
	StringC directory;
	SArray1dC<Point2dC> a_mean;
	Affine2dC sim_proj;
	SArray1dC<Point2dC> ProjectAffinePoints(SArray1dC<Point2dC> const &pts, Affine2dC const &proj);
};

#endif
