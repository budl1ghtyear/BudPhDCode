#ifndef CHIHISTMEAS_HH
#define CHIHISTMEAS_HH
#include "Ravl/SArray1dIter2.hh"
#include "Ravl/Array2dIter.hh"
#include "BaseMeasurementModel.hh"
#include "TrackingUtilityFunctions.hh"
#include "UtilityFunctions.hh"
#include "Ravl/PatternRec/Sample.hh"
#include "Ravl/PatternRec/SampleIter.hh"
#include "Ravl/LinePP2d.hh"
#include "Ravl/RealRange2d.hh"
#include "Ravl/Polygon2d.hh"
#include "Ravl/RealHistogram3d.hh"
#include "Ravl/Point3d.hh"
#include "opencv/highgui.h"
#include "opencv/cv.h"
using namespace cv;
using namespace RavlN;
using namespace RavlImageN;

class ChiHistMeasurementModelC : public BaseMeasurementModelC
{
	public:
		ChiHistMeasurementModelC(){}
		~ChiHistMeasurementModelC(){}
		//Construct the measurement model using a simple mouth-region image from which we can obtain edge-measurements
		ChiHistMeasurementModelC(ImageC<RealRGBValueC> const &img, ImageRectangleC const &rect): im(img.Copy()),imrect(rect)
		{}
		RealT Measure(ParticleC &pt);
		RealT MeasureQuality(SArray1dC<LinePP2dC> &norms);
		RealT ChiDistance(ImageRectangleC const &rec1,ImageRectangleC const &rec2, Point3dC const &mask = Point3dC(8.,8.,8.));
		//Accessor Methods:
		ImageC<RealRGBValueC> GetOriginalImage(void) const  {return im;}
		ImageRectangleC GetPrevRect(void) const {return imrect;}
		void SetPrevRect(ImageRectangleC const &rect) {imrect = rect;}
		ImageRectangleC GetImageRectangleFromPts(SArray1dC<Point2dC> const &pts);
	protected:
		ImageC<RealRGBValueC> im; //original colour image
		ImageRectangleC imrect;
		template <typename T>
		//void GetHistogram(ImageC<T> const &img, ImageRectangleC const &rect, CvHistogram* hist, Point3dC const &mask = Point3dC(8.,8.,8.))
		//~ RealHistogram3dC GetHistogram(ImageC<T> const &img, ImageRectangleC const &rect, Point3dC const &mask = Point3dC(8.,8.,8.))
		SArray3dC<UIntT> GetHistogram(ImageC<T> const &img, ImageRectangleC const &rect, Point3dC const &mask = Point3dC(8.,8.,8.))
		{
			ImageC<T> myimg(img,rect);
			//~ Index3dC step(256./mask.X(),256./mask.Y(),256./mask.Z());
			//~ RealHistogram3dC myhist(Point3dC(0.,0.,0.),Point3dC(mask.X(),mask.Y(),mask.Z()),step);
			// make sure that the histogram has proper size and type
			//int histSize[] = {mask.X(),mask.Y(),mask.Z()};
			//hist.create(3, histSize, CV_32F);
			// and clear it
			//hist = Scalar(0);
			//int histSize[] = {mask.X(),mask.Y(),mask.Z()};
			//float x_ranges[] = { 0, mask.X() };
			//float y_ranges[] = { 0, mask.Y() };
			//float z_ranges[] = { 0, mask.Z() };
			//float* ranges[] = {x_ranges,y_ranges,z_ranges};
			//hist = cvCreateHist( 3, histSize, CV_HIST_ARRAY, ranges, 1 );			
			SArray3dC<UIntT> myhist(mask.X(),mask.Y(),mask.Z());
			for(Array2dIterC<T>it(myimg); it; it++)
			{
				T pix = *it;
				// we could have incremented the cells by 1.f/(image.rows*image.cols)
				// instead of 1.f to make the histogram normalized.
				//hist.at<float>(pix[0]*mask.X()/256, pix[1]*mask.Y()/256, pix[2]*mask.Z()/256) += 1.f;
				//*cvGetHistValue_3D(hist,pix[0]*mask.X()/256, pix[1]*mask.Y()/256, pix[2]*mask.Z()/256) += 1.f;
				Point3dC pt(pix[0],pix[1],pix[2]);
				//myhist.Vote(pt,1.);
				myhist[Floor(pix[0]*mask.X()/256.)][Floor(pix[1]*mask.Y()/256.)][Floor(pix[2]*mask.Z()/256.)] += 1;
			}
			return myhist;
		}
};

#endif
