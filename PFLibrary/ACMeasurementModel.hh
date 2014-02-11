#ifndef ACMeasurementModelC_HH
#define ACMeasurementModelC_HH

#include "BaseMeasurementModel.hh"
#include "Ravl/SArray1d.hh"
#include "Ravl/SArray2dIter.hh"
#include "Ravl/Array1d.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/SArray2dIter2.hh"
#include "Ravl/SArray2dIter3.hh"
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
#include "BasisSpline.hh"
#include "UtilityFunctions.hh"
//This class implements the Active Contours Based measurement model as detailed by Montse Pardas in her paper:
//Motion estimation based tracking of active contours
//It includes the estimation of an internal and external energy term as described inside the paper

using namespace RavlN;
using namespace RavlImageN;
class ACMeasurementModelC : public BaseMeasurementModelC
{
	public:
	ACMeasurementModelC(){}
	virtual ~ACMeasurementModelC(){}
	ACMeasurementModelC(const ImageC<RealT> &like_vec_img, const ImageC<Vector2dC> &img, const SArray1dC<Point2dC> &prev, const RealT &g = 0.5, const RealT &l= 0.5, const RealT &b= 1, const UIntT &w= 8):vec_img(like_vec_img.Copy()),oflow_img(img.Copy()),prev_snaxels(prev.Copy()),gamma(g),lamda(l),beta(b),window(w) {}
	RealT Measure(ParticleC &pt);
	protected:
	RealT ComputeInternalEnergy(const SArray1dC<Point2dC> &cur) const ;
	RealT ComputeExternalEnergy(const SArray1dC<Point2dC> &cur, const BasisSplineC &bspl);
	RealT ComputeContourEnergy(const SArray1dC<Point2dC> &cur, const BasisSplineC &bspl);
	RealT ComputeMotionCompensationTerm(const SArray1dC<Point2dC> &cur) ;
	ImageRectangleC RangeCheck(ImageRectangleC &ind1, ImageRectangleC &ind2)
	{
		//ind1 is the proposed rectangle and ind2 is the image rectangle within which the ind1 should be
		ImageRectangleC rect;
		//if the ind1.TRow() is positionally lower than the ind2.TRow(), return it
		rect.TRow() = (ind2.TRow() < ind1.TRow() ? ind1.TRow() : ind2.TRow());
		//if the ind1.LCol() is positionally more to the left than the ind2.LCol(), return it
		rect.LCol() = (ind2.LCol() < ind1.LCol() ? ind1.LCol() : ind2.LCol());
		//if the ind1.BRow() is positionally higher than the ind2.BRow(), return it
		rect.BRow() = (ind2.BRow() < ind1.BRow() ? ind2.BRow() : ind1.BRow());
		//if the ind1.RCol() is positionally more to the right than the ind2.RCol(), return it
		rect.RCol() = (ind2.RCol() < ind1.RCol() ? ind2.RCol() : ind1.RCol());
		rect.TRow() = rect.TRow() < rect.BRow() ? rect.TRow() : rect.BRow();
		rect.LCol() = rect.LCol() < rect.RCol() ? rect.LCol() : rect.RCol();
		return rect;
	}
	ImageC<RealT> vec_img;
	ImageC<Vector2dC> oflow_img;
	SArray1dC<Point2dC> prev_snaxels;
	RealT gamma, lamda, beta; //Heuristic constants to be used inside the equations
	UIntT window;
};

#endif
