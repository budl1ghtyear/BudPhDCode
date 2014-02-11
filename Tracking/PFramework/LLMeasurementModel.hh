#ifndef LLMEASUREMENTMODEL_HH
#define LLMEASUREMENTMODEL_HH

////////////////////////////////
// Name - JMeasurementModelC
// Author - Bud Goswami
// Notes : Implements the log-likelihood measurement model based on URS thesis Ch 4 - using J function
////////////////////////////////

#include "BaseMeasurementModel.hh"
#include "Ravl/Array1d.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/Array2dIter2.hh"
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
using namespace RavlN;
using namespace RavlImageN;
using namespace RavlConstN;

class LLMeasurementModelC : public BaseMeasurementModelC
{
public:
	LLMeasurementModelC() {}
	LLMeasurementModelC(const ImageC<RealT> &jim, const MeanCovarianceC &lp, const MeanCovarianceC &sk, const UIntT &win);
	LLMeasurementModelC(const ImageC<RealT> &jim, const MeanCovarianceC &lp, const MeanCovarianceC &sk);
	RealT Measure(ParticleC &pt);
	
protected:
	RealT ComputeWeight(const Array1dC<Point2dC> &ar);
	RealT ComputeWeight(const Array1dC<Point2dC> &ar, const BSplineC &bspl);
	RealT doComputeWeight(const Array1dC<Point2dC> &ar, const BSplineC &bspl);
	//RealT GetGaussianValue(const RealT &sigma, const RealT &xval);
	ImageC<RealT> vec_img;
	MeanCovarianceC lip;
	MeanCovarianceC skin;
	UIntT window;
};

#endif
