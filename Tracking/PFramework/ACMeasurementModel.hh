#ifndef ACMEASUREMENTMODEL_HH
#define ACMEASUREMENTMODEL_HH

#include "BaseMeasurementModel.hh"
#include "Ravl/Array1d.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/Array2dIter2.hh"
#include "Ravl/Array2dIter3.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Image/EdgeSobel.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Index2d.hh"
#include "Ravl/LinePP2d.hh"
#include "Ravl/Math.hh"
#include "Ravl/MeanCovariance.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/Vector.hh"
#include "BSplineC.hh"

//This class implements the Active Contours Based measurement model as detailed by Montse Pardas in her paper:
//Motion estimation based tracking of active contours
//It includes the estimation of an internal and external energy term as described inside the paper

using namespace RavlN;
class ACMeasurementModel : public BaseMeasurementModelC
{
	ACMeasurementModel(){}
	virtual ~ACMeasurementModel(){}
	ACMeasurementModel(const ImageC<RealT> &like_vec_img, const ImageC<Vector2dC> &img, const Array1dC<Point2dC> &prev, const RealT &g, const RealT &l, const RealT &b, const UIntT &w);
	RealT Measure(ParticleC &pt);
	protected:
	RealT ComputeInternalEnergy(const Array1dC<Point2dC> &cur) const ;
	RealT ComputeExternalEnergy(const Array1dC<Point2dC> &cur, const BSplineC &bspl) const ;
	RealT ComputeContourEnergy(const Array1dC<Point2dC> &cur, const BSplineC &bspl) const;
	RealT ComputeMotionCompensationTerm(const Array1dC<Point2dC> &cur) const;
	
	ImageC<RealT> vec_img;
	ImageC<Vector2dC> oflow_img;
	Array1dC<Point2dC> prev_snaxels;
	RealT gamma, lamda, beta; //Heuristic constants to be used inside the equations
	UIntT window;
}

#endif
