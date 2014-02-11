#ifndef PCAPROJECTION_HH
#define PCAPROJECTION_HH
//Class to perform PCA Projection both to and from PCA space
#include "Ravl/Array2d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/SArray1d.hh"
#include "Ravl/SArray1dIter.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/Matrix.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/PatternRec/FuncMeanProjection.hh"
#include "Ravl/PatternRec/FuncLinear.hh"
#include "Ravl/Vector.hh"
#include "UtilityFunctions.hh"
using namespace RavlN;
using namespace RavlImageN;

class PCAProjectionC
{
	public:
	PCAProjectionC(){}
	~PCAProjectionC(){}
	
	PCAProjectionC(StringC const &pcadir, UIntT const &numcp=10); //: Constructor to use with a specified input directory path
	PCAProjectionC(VectorC const &mean, MatrixC const &eigenvect, VectorC const &eigenval, UIntT const &numcp=10); //:Standard constructor to use
	PCAProjectionC(VectorC const &mean, MatrixC const &eigenvect, VectorC const &eigenval, DirectoryC const &dir, UIntT const &numcp=10):pca_mean(mean),pca_loadings(eigenvect),pca_sdev(eigenval),directory(dir),numcomp(numcp) {;}
	PCAProjectionC(PCAProjectionC const &proj)
	{
		this->pca_mean = proj.GetMean();
		this->pca_loadings = proj.GetEigenVectors();
		this->pca_sdev = proj.GetEigenValues();
		this->directory = proj.GetDirectory();
		this->numcomp = proj.GetNumComponents();
		//cout<<"PCAProjectionC Copy Called"<<endl;
	}
	PCAProjectionC Copy() const
	{
		return PCAProjectionC(pca_mean,pca_loadings,pca_sdev,directory,numcomp);
	}
	//Useful Methods
	VectorC PerformPCAProjection(SArray1dC<Point2dC> const &points);
	VectorC PerformPCAProjection(Array1dC<Point2dC> const &points);
	SArray1dC<Point2dC> PerformInverseProjection(VectorC const &vect); 
	VectorC GetOuterBound(UIntT const &n=1)
	{
		//PCA_SDEV is the eigenvalue matrix obtained from R
		//It contains the roots of the individual eigenvalues and therefore is an estimate of the standard deviation of the projection anyway.
		//We perform window based perturbation inside the propagation model
		VectorC eigs = pca_sdev.From(0,numcomp).Copy();
		VectorC res = eigs * n;
		//cout<<"GetOuterBound() in PCAProjectionC returns - "<<res<<endl;
		return res;
		//return (pca_sdev*n);
	}
	//Accessor methods:
	VectorC GetMean(void) const {return pca_mean;}
	VectorC GetEigenValues(void) const {return pca_sdev;}
	MatrixC GetEigenVectors(void) const {return pca_loadings;}
	DirectoryC GetDirectory(void) const {return directory;}
	UIntT GetNumComponents(void) const {return numcomp;}	
	protected:
	VectorC pca_mean;
	MatrixC pca_loadings;
	VectorC pca_sdev;
	StringC directory;
	UIntT numcomp;
};



#endif
