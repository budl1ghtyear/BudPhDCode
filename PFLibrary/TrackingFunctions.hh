//TrackingFunctions.hh
//Author - Bud Goswami
//Purpose - to create a header file that contains the function definitions of mouth-region tracking
//Use RAVL as interface to OpenCV face detection and OMNI Eye Centre Locator
#ifndef TrackingFunctions_HH
#define TrackingFunctions_HH
#include "Ravl/Array1d.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/Affine2d.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Ellipse2d.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/DrawFrame.hh"
#include "Ravl/Image/DrawLine.hh"
#include "Ravl/Image/OpenCVConvert.hh"
#include "Ravl/Image/ImageRectangle.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/IO.hh"
#include "Ravl/MeanCovariance.hh"
#include "Ravl/Option.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/Tuple2.hh"
#include "Ravl/SArray1d.hh"
#include "Ravl/SArray1dIter.hh"
#include "Ravl/String.hh"
#include "Ravl/SumsNd2.hh"
#include "Ravl/Vector.hh"
#include "Ravl/Image/GaussConvolve2d.hh"
#include "Ravl/Image/WarpAffine.hh"
#include "Ravl/Image/HistogramEqualise.hh"
#include "BSplineC.hh"
using namespace RavlImageN;
using namespace RavlConstN;

Tuple2C<ImageC<RealT>,Affine2dC> ConvertImages(const ImageC<RealT> & image, const Point2dC &nle,const Point2dC &nre);
Array1dC<Point2dC> ConvertVecToPointArray(const VectorC &v);
ImageC<RealT> ConvertRGBToRealImage(const ImageC<RealRGBValueC> &im);
Array1dC<Index2dC> GetEpsBoundary(const Array1dC<Index2dC> &ar, ImageC<RealT> &mc, ImageC<RealT> &sk);
Array1dC<Point2dC> GetLipContourPoints(const DListC<Index2dC> &lp, const UIntT &res, const ImageC<UIntT> &img, const RealT &search, const ImageC<RealT> &jim);
//TEMPLATED FUNCTIONS
template <typename T>
DListC<Index2dC> ClusterIdent(DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > &finalpop, const ImageC<T> &img, const IntT &segtype)
{
	DListC<Index2dC> lip;
	//Method just chooses the type of function to use for cluster identification depending on segtype
	switch(segtype)
	{
		case 0:
		{
			lip = DensityCG(finalpop, img);
			break;
		}
		case 1:
		{
			lip = JCC(finalpop, img);
			break;
		}
	}
	return lip;
}
template <typename T>
MeanCovarianceC GetMC(const DListC<Index2dC> &pts, const ImageC<T> &img)
{
	UIntT dim = img[img.Frame().TopLeft()].Size();
	SumsNd2C sum(dim);
	for(DLIterC<Index2dC> d(pts); d; d++)
	{
		VectorC vec(dim);
		T type = img[(*d)];
		for(UIntT i = 0; i < dim; i++)
		{
			vec[i] = type[i];
		}
		sum += vec.Copy();
	}
	return sum.MeanCovariance();
}
template <typename T>
ImageC<RealT> GetMSImg(const ImageC<T> &img, const MeanCovarianceC &mc, const MeanCovarianceC &sk)
{
	UIntT dim = img[img.Frame().TopLeft()].Size();
	ImageC<RealT> out(img.Frame());
	out.Fill(0.0);
	for(Array2dIter2C<T, RealT> it(img,out); it; it++)
	{
		VectorC vec(dim);
		for(UIntT i = 0; i < dim; i++)
		{
			vec[i] = it.Data1()[i];
		}
		
		RealT val = mc.MahalanobisDistance(vec);
		if(val < 6.0)
			it.Data2() = val;
		else if ((val < sk.MahalanobisDistance(vec))&&(val>= 6.0))
			it.Data2() = 6.0;
		else
			it.Data2() = 0.0;
	}
	return out;
}
template <typename T>
ImageC<RealT> ComputeJ(const ImageC<T> &img, const MeanCovarianceC &sk, const MeanCovarianceC &lp)
{
	RealT log_sk = Log(sk.Covariance().Det()); RealT log_lp = Log(lp.Covariance().Det());  
	UIntT size = img[img.Frame().TopLeft()].Size();
	ImageC<RealT> res(img.Frame(),0.0);
	for(Array2dIter2C<T,RealT> it(img,res); it; it++)
	{
		VectorC vec(size);
		for(UIntT i = 0; i < size; i++)
		{
			vec[i] = it.Data1()[i];
		}
		RealT j_sk = -log_sk - sk.MahalanobisDistance(vec);
		RealT j_lp = -log_lp - lp.MahalanobisDistance(vec);
		it.Data2() = j_lp - j_sk;
	}
	return res;
}
//TESTING FUNCTIONS
void ComputeTest(const SArray1dC<Point2dC> &tres,BSplineC &spl, const ImageC<RealRGBValueC> &rgbimg, const UIntT &res);
#endif

