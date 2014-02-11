#ifndef LIPCLUSTERC_HH
#define LIPCLUSTERC_HH
//////////////////////////////////////////////////////////////////////////////
//! file =      "Ravl/LipCode/Segmentation/SegmentationLibrary/LipClustering.hh"   
//! author =    "Bud Goswami"
//! lib =       SegLib
//! date =      "14.08.09"
//! userlevel = Normal
//! docentry =  "Ravl.API.LipCode.Segmentation.SegmentationLibrary"


//INCLUDE FILES
#include "Ravl/Array1d.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/ConnectedComponents.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/DrawFrame.hh"
#include "Ravl/Image/OpenCVConvert.hh"
#include "Ravl/Image/ImageRectangle.hh"
#include "Ravl/Image/Segmentation.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/PatternRec/DesignGaussianMixture.hh"
#include "Ravl/PatternRec/Sample.hh"
#include "Ravl/PatternRec/DistanceMahalanobis.hh"
#include "Ravl/PatternRec/Classifier.hh"
#include "Ravl/PatternRec/DesignKMeans.hh"
//#include "Ravl/PatternRec/DesignFuzzyCMeansCluster.hh"
#include "DesignFuzzyCMeansCluster.hh"
//BUDSW
#include "CascadedMCD.hh"
#include "KSmirnovMCD.hh"
#include "MouthRegion.hh"
using namespace RavlN;
using namespace RavlImageN;
using namespace OmniN;

typedef	enum {CMCD=0,KSMCD=1,KM=2,FCM=3}ClusteringType; 
//:Enumerated type to specify the sort of clustering you want to perform on an image<p>
//:for CascadedMCD, Kolmorogorov-Smirnov MCD, K-Means and Fuzzy C Means<p>
//: This is global so that the calling program can directly using it without any Integer to Enumerated type conversions


//:Class encapsulates the clustering mechanism from a user. Use this class as the calling program for all your image clustering requirements
template <class T>
class LipClusteringC
{
	public:
	LipClusteringC(){};
	~LipClusteringC(){};
	//My constructors
	//LipClusteringC(const ImageC<T> &img, const RealT &hval = (RealT)0.75, const UIntT &numclust = (UIntT) 3);
	LipClusteringC(const ImageC<T> &img, const RealT &hval, const UIntT &numclust);
	//:Constructor with the image, a h confidence measure for MCD based clustering and number of clusters for KMeans and FCM based algorithms
	//: Constructor with some image type(NOTE: Image needs to be multidimensional), h value for MCD and number of clusters for KMeans or Fuzzy C Means 
	DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > Apply(const ClusteringType &ctype);	
	//: Apply() method with clustering type option (look at enumerated constants) to perform appropriate clustering operation
	
	//Accessor Methods
	ImageC<T> GetOriginalImage(void) const {return inp_img;}
	//:GetOriginalImage()
	RealT GetHValue(void) const {return h;}
	//:GetHValue()
	UIntT GetNumClusters(void) const {return n;}
	//:GetNumClusters()
	protected:
	ImageC<T> inp_img;
	RealT h; //:hvalue for MCD algorithms (default = 0.75)
	UIntT n; //:number of clusters for KM and FCM (default = 3)
};

//#include "LipClusteringC.cc"

template <class T>
LipClusteringC<T>::LipClusteringC(const ImageC<T> &img, const RealT &hval = (RealT)0.75, const UIntT &numclust = (UIntT) 3):h(hval),n(numclust)
{
	inp_img = img.Copy();
}


template <class T>
DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > LipClusteringC<T>::Apply(const ClusteringType &ctype)
{
	DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > result;
	switch(ctype)
	{
		case CMCD:
		{
			CascadedMCD mcd(inp_img,h);
			//Get the data we want to cluster and input into the clustering algorithm
			result = mcd.Apply();
			break;
		}
		case KSMCD:
		{
			KSmirnovMCD mcd(inp_img,h); //load it up with the image type we would like
			//Get the data we want to cluster and input into the clustering algorithm
			result = mcd.Apply();
			break;	
		}
		case KM:
		{	
			DesignKMeansC km(n);
			SampleC<VectorC> smpl = GetSample(inp_img);
			ClassifierC cl = km.Apply(smpl);
			ImageC<UIntT> clustimg(inp_img.Frame(),0);
			UIntT dim = inp_img[inp_img.Frame().TopLeft()].Size();
			DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > out;
			SArray1dC<DListC<Tuple2C<VectorC,Index2dC> > > pix(n);
			for(Array2dIterC<T> it(inp_img); it; it++)
			{
				VectorC v(dim);
				for(UIntT i = 0; i < dim; i++)
				{
					v[i] = (*it)[i];
				}
				UIntT label = cl.Classify(v);
				clustimg[it.Index()] = label;
				//Insert VectorC and Index into the output cluster
				Tuple2C<VectorC,Index2dC> tpl(v.Copy(), it.Index());
				pix[label].Append(tpl);
			}
			//Now obtain the MeanCovariance Values for each cluster
			for(SArray1dIterC<DListC<Tuple2C<VectorC,Index2dC> > > it(pix); it; it++)
			{
				SumsNd2C sum(dim);
				bool SampleStats = false;
				for(DLIterC<Tuple2C<VectorC,Index2dC> > it2((*it)); it2; it2++)
				{
					sum += (*it2).Data1();
				}
				Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > elem(sum.MeanCovariance(SampleStats).Copy(),(*it));
				out.Append(elem.Copy());
			}
			result = out.Copy();	
			break;		
		}
		case FCM:
		{
			MyDesignFuzzyCMeansClusterC km(n);
			SampleC<VectorC> smpl = GetSample(inp_img);
			ClassifierC cl = km.Apply(smpl);
			ImageC<UIntT> clustimg(inp_img.Frame(),0);
			UIntT dim = inp_img[inp_img.Frame().TopLeft()].Size();
			DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > out;
			SArray1dC<DListC<Tuple2C<VectorC,Index2dC> > > pix(n);
			for(Array2dIterC<T> it(inp_img); it; it++)
			{
				VectorC v(dim);
				for(UIntT i = 0; i < dim; i++)
				{
					v[i] = (*it)[i];
				}
				UIntT label = cl.Classify(v);
				clustimg[it.Index()] = label;
				//Insert VectorC and Index into the output cluster
				Tuple2C<VectorC,Index2dC> tpl(v.Copy(), it.Index());
				pix[label].Append(tpl);
			}
			//Now obtain the MeanCovariance Values for each cluster
			for(SArray1dIterC<DListC<Tuple2C<VectorC,Index2dC> > > it(pix); it; it++)
			{
				SumsNd2C sum(dim);
				bool SampleStats = false;
				for(DLIterC<Tuple2C<VectorC,Index2dC> > it2((*it)); it2; it2++)
				{
					sum += (*it2).Data1();
				}
				Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > elem(sum.MeanCovariance(SampleStats).Copy(),(*it));
				out.Append(elem.Copy());
			}
			result = out.Copy();
			break;					
		}
		default:
		{
			result = Apply(CMCD) ;
			break;
		}
	}
	return result;	
}

#endif
