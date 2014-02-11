#ifndef EDGEMEASUREMENTMODEL_HH
#define EDGEMEASUREMENTMODEL_HH

#include "BaseMeasurementModel.hh"
#include "TrackingUtilityFunctions.hh"
#include "UtilityFunctions.hh"
#include "Ravl/PatternRec/Sample.hh"
#include "Ravl/PatternRec/SampleIter.hh"
#include "Ravl/LinePP2d.hh"
using namespace RavlN;
using namespace RavlImageN;

class EdgeMeasurementModelC : public BaseMeasurementModelC
{
	public:
	EdgeMeasurementModelC(){}
	//Construct the measurement model using a simple mouth-region image from which we can obtain edge-measurements
	EdgeMeasurementModelC(ImageC<RealRGBValueC> const &img, ImageRectangleC const &imrec, IntT const &win = 10);
	RealT Measure(ParticleC &pt);
	RealT MeasureQuality(SArray1dC<LinePP2dC> &curvenorms);
	//Accessor Methods:
	ImageC<RealRGBValueC> GetOriginalImage(void) const  {return im;}
	ImageC<RealT> GetEdgeImage(void) const  {return edge;}
	UIntT GetSearchWindow(void) const  {return searchwindow;}
	protected:
	ImageC<RealRGBValueC> im; //original colour image
	ImageC<RealT> edge; //edge image
	ImageRectangleC myimgrect;
	IntT searchwindow; //length of edge search window
};

#endif
