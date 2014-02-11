#ifndef EDGEMEASUREMENTMODEL_HH
#define EDGEMEASUREMENTMODEL_HH

#include "BaseMeasurementModel.hh"
#include "TrackingUtilityFunctions.hh"
#include "UtilityFunctions.hh"
using namespace RavlN;
using namespace RavlImageN;

class EdgeMeasurementModelC : public BaseMeasurementModelC
{
	public:
	EdgeMeasurementModelC(){}
	//Construct the measurement model using a simple mouth-region image from which we can obtain edge-measurements
	EdgeMeasurementModelC(ImageC<RealRGBValueC> const &img, IntT const &win = 10);
	RealT Measure(ParticleC &pt);
	//Accessor Methods:
	ImageC<RealRGBValueC> GetOriginalImage(void) const  {return im;}
	ImageC<RealT> GetEdgeImage(void) const  {return edge;}
	UIntT GetSearchWindow(void) const  {return searchwindow;}
	protected:
	ImageC<RealRGBValueC> im; //original colour image
	ImageC<RealT> edge; //edge image
	IntT searchwindow; //length of edge search window
};

#endif
