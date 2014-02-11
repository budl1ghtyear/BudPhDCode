#ifndef TrackingUtilityFunctions_HH
#define TrackingUtilityFunctions_HH

//Class to perform the utility functions required for tracking
#include "Ravl/Affine2d.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "BasisSpline.hh"
#include "Particle.hh"
#include "SystematicResampling.hh"
#include "SIR.hh"
#include "Ravl/Image/LMSOpticFlow.hh"
using namespace RavlN;
using namespace RavlImageN;

//Tuple2C<BasisSplineC,Affine2dC> GetLipPoints(FilenameC const &lipfile);
//Tuple2C<BasisSplineC,Affine2dC> GetLipPoints(ImageC<RealRGBValueC> const &lipfile);
//ImageC<UIntT> MakeFullSizeImage(const ImageC<UIntT> &r_img, const ImageRectangleC &frame);
//Tuple2C<ImageC<RealT>, ImageC<Vector2dC> > GetMeasurementImages(ImageC<RealRGBValueC> const &im1, Tuple2C<Point2dC,Point2dC> const &eyes1,ImageC<RealRGBValueC> const &im2, Tuple2C<Point2dC,Point2dC> const &eyes2);
ImageC<RealT> GetSusanImage(ImageC<RealT> const &img);
ImageC<RealT> ConvertRGBToRealImage(const ImageC<RealRGBValueC> &im);
//SArray1dC<Point2dC> AffineTransformPoints(const SArray1dC<Point2dC> &pts, const Affine2dC &aff);
void DisplayImage(const SArray1dC<Point2dC> &pts, const UIntT &size=2);
ImageRectangleC ObtainMouthLocalisationWindow(SArray1dC<Point2dC> const &pts, UIntT const &src);
template <typename T>
ImageC<T> NormaliseImage(ImageC<T> const &im)
{
	ImageC<T> result(im.Frame());
	T maxval = im[im.IndexOfMax()];
	T minval = im[im.IndexOfMin()];
	T den = maxval - minval;
	for(Array2dIter2C<T,T> it(im, result); it; it++)
	{
		it.Data2() = (it.Data1() - minval)/den;
	}	
	return result;
}
#endif
