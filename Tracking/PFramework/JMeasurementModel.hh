#ifndef JMEASUREMENTMODEL_HH
#define JMEASUREMENTMODEL_HH

////////////////////////////////
// Name - JMeasurementModelC
// Author - Bud Goswami
// Notes : Implements the measurement model based on URS thesis Ch 4 - using J function
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

class JMeasurementModelC : public BaseMeasurementModelC
{
public:
	JMeasurementModelC(){}
	~JMeasurementModelC(){}
	JMeasurementModelC(const MeanCovarianceC &lp, const MeanCovarianceC &sk, const ImageC<RealRGBValueC> &img, const IntT &wd);
	virtual RealT Measure(ParticleC &pt);
	virtual void operator=(const JMeasurementModelC &oth)
	{
		lip = oth.GetLipMC();
		skin = oth.GetSkinMC();
		origimg = oth.GetImageSpace();
		window = oth.GetWindow();
	}
	//Accessor Methods
	MeanCovarianceC GetLipMC(void) const {return lip;}
	MeanCovarianceC GetSkinMC(void) const {return skin;}
	ImageC<RealRGBValueC> GetImageSpace(void) const {return origimg;}
	IntT GetWindow(void) const {return window;}
protected:
	ImageC<VectorC> GetNormRGVals(const ImageC<RealRGBValueC> &im);
	Array1dC<Index2dC> VectortoPoint(const VectorC &vec);
	Array1dC<Point2dC> VectortoPointArray(const VectorC &vec);
	ImageC<RealT> ComputeJImage(const ImageC<RealRGBValueC> &img, const MeanCovarianceC &mc);
	ImageC<RealT> ComputeJImage(const ImageC<VectorC> &img, const MeanCovarianceC &mc);
	RealT GetJWeight(const Array1dC<Point2dC> &ar, ImageC<RealT> &mc, ImageC<RealT> &sk);
	RealT GetJWeight(const Array1dC<Point2dC> &ar);
	//RealT GetGaussianValue(const RealT &sigma, const RealT &xval);
	RealT ComputeJ(const MeanCovarianceC &mc, const VectorC &vec);
	MeanCovarianceC lip;
	MeanCovarianceC skin;
	ImageC<RealRGBValueC> origimg;
	IntT window;
};

#endif
