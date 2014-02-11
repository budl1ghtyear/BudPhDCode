#ifndef ICAPROJECTION_HH
#define ICAPROJECTION_HH
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

class ICAProjectionC
{
	public:
	ICAProjectionC(){}
	
	ICAProjectionC(DirectoryC const &dir);
	ICAProjectionC(MatrixC const &w, MatrixC const &a):mixing(a),unmixing(w){;}
	ICAProjectionC(ICAProjectionC const &ica)
	{
		this->mixing = ica.GetMixingMatrix();
		this->unmixing = ica.GetUnmixingMatrix();
	}
	ICAProjectionC Copy() const 
	{
		return ICAProjectionC(unmixing,mixing);
	}
	//Useful Methods:
	VectorC PerformICAProjection(SArray1dC<Point2dC> const &points);
	VectorC PerformICAProjection(Array1dC<Point2dC> const &points);
	SArray1dC<Point2dC> PerformInverseProjection(VectorC const &vect); 
	//Accessor Methods:
	MatrixC GetMixingMatrix(void) const {return mixing;}
	MatrixC GetUnmixingMatrix(void) const {return unmixing;}
	VectorC GetOuterBound(UIntT const &n=1)
	{
		VectorC eigs(unmixing.Cols());
		eigs.Fill(n);
		return eigs;
		//return (pca_sdev*n);
	}		
	protected:
	MatrixC mixing;
	MatrixC unmixing;
};


#endif
