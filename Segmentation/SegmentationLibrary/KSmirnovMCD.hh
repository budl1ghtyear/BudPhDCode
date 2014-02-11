/*
File: KSmirnovMCD.hh
 
Author: Bud Goswami
Date: 16.02.09
*/

#ifndef KSmirnovMCD_HH
#define KSmirnovMCD_HH
//////////////////////////////////////////////////////////////////////////////
//! file =      "Ravl/LipCode/Segmentation/SegmentationLibrary/KSmirnovMCD.hh"   
//! author =    "Bud Goswami"
//! lib =       SegLib
//! date =      "16.02.09"
//! userlevel = Normal
//! docentry =  "Ravl.API.LipCode.Segmentation.SegmentationLibrary"
//INCLUDE FILES
#include "Ravl/Array1dIter2.hh"
#include "Ravl/SumsNd2.hh"

//INCLUDE OWN SW
#include "chisquaredistr.hh"
#include "MCD.hh"


using namespace RavlN;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------
// **********  KSMIRNOVMCD  **********************************************
// ---------------------------------------------------------------------------
//:Class - declares the interface to the KSmirnov Minimum Covariance Determinant Function<p>
//:Takes in as an input, a set of Sample of Vectors. <p>
//:Performs iterative MCD estimation on that set of vectors until it runs out of data<p>
//:Can return various final estimates from procedure.<p>
class KSmirnovMCD : public MCD
{
	//Public Methods:
	//Default Constructor and Destructor
public:
	KSmirnovMCD() { cout<<"Empty KSmirnovMCD Object called"<<endl;}
	~KSmirnovMCD() {}
	
	
	KSmirnovMCD(const ImageC<RealRGBValueC> &rgbimg, const RealT &hval, const UIntT &numiter) : MCD(rgbimg,hval,numiter), minsamp(100),maxclust(4){}
	//:Constructor inputs: ImageC<RealRGBValueC>, RealT (containing h) and UIntT (containing the number of iterations)<p>
	//:variations of above (Constructor Must Always Have At Least ImageC<RealRGBValueC> object)<p>
	//:More image variations will require recoding<p>
	//:DEFAULTS VALUES:	<p>
	//:h = 0.75 numiter = 10 minsamp =100 <p>
	KSmirnovMCD(const ImageC<RealRGBValueC> &rgbimg, const RealT &hval):MCD(rgbimg, hval), minsamp(100),maxclust(4){}
	KSmirnovMCD(const ImageC<RealRGBValueC> &rgbimg, const UIntT &numiter): MCD(rgbimg, numiter), minsamp(100),maxclust(4){}
	KSmirnovMCD(const ImageC<RealRGBValueC> &rgbimg) : MCD(rgbimg), minsamp(100),maxclust(4){}
	
	KSmirnovMCD(const ImageC<VectorC> &vecimg, const RealT &hval, const UIntT &numiter) : MCD(vecimg, hval, numiter), minsamp(100),maxclust(4){}
	//:In case you want an image with just an arbitrary vector type
	KSmirnovMCD(const ImageC<VectorC> &vecimg, const RealT &hval) : MCD(vecimg,hval), minsamp(100),maxclust(4){}
	KSmirnovMCD(const ImageC<VectorC> &vecimg, const UIntT &numiter): MCD(vecimg, numiter), minsamp(100),maxclust(4){}
	KSmirnovMCD(const ImageC<VectorC> &vecimg) : MCD(vecimg), minsamp(100),maxclust(4){}	
	
	KSmirnovMCD(SampleData &smpl, const RealT &hval, const UIntT &numiter):MCD(smpl,hval), minsamp(100),maxclust(4){}
	//:Constructors with SampleC<VectorC>, RealT and UIntT
	//:NOTE that using SampleDataC for images assumes that the vector includes index information
	//:Form: vec[0] = imgdata0...vec[dim-1] = imgdata_dim; vec[dim] = Index.Row(), vec[dim+1] = Index.Col() 
	KSmirnovMCD(SampleData &smpl, const RealT &hval):MCD(smpl, hval), minsamp(100),maxclust(4){}
	KSmirnovMCD(SampleData &smpl, const UIntT &numiter) : MCD(smpl, numiter), minsamp(100),maxclust(4){}
	KSmirnovMCD(SampleData &smpl) : MCD(smpl), minsamp(100),maxclust(4){}
	
	KSmirnovMCD(const ImageC<TFVectorC<RealT,3> > &rgbimg, const RealT &hval, const UIntT &numiter):MCD(rgbimg,hval,numiter), minsamp(100) {}
	//:For usage with other 3D Real Colour Spaces 
	KSmirnovMCD(const ImageC<TFVectorC<RealT,3> > &rgbimg, const RealT &hval):MCD(rgbimg,hval), minsamp(100) {}
	KSmirnovMCD(const ImageC<TFVectorC<RealT,3> > &rgbimg, const UIntT &numiter):MCD(rgbimg, numiter), minsamp(100) {}
	KSmirnovMCD(const ImageC<TFVectorC<RealT,3> > &rgbimg):MCD(rgbimg), minsamp(100) {}	
	//Accessor Functions
	void SetStartPop(const DListC<Tuple2C<VectorC,Index2dC> > &pop) 
	{
		start_pop.Empty();
		start_pop = pop.Copy();
	}
	void SetMinSamp(const UIntT &ms) 
	{
		if (ms >= 0)
			minsamp = ms;
		else
			cerr<<"MIN SAMPLES must be greater than zero";
	}
	UIntT GetMinSamp(void) const {return minsamp;}
	void SetMaxClust(const UIntT &ms) 
	{
		if (ms >= 0)
			maxclust = ms;
		else
			cerr<<"MAX CLUSTERS must be greater than zero";
	}
	UIntT GetMaxClusters(void) const {return maxclust;}
	
	DListC< Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > Apply(void);
	//:Apply Method
protected:
	//Private Data Members
	UIntT minsamp; //: Store the minimum number of samples, default value is 100
	UIntT maxclust; //: Variable to stop the algorithm once it extracts maxclust number of clusters.
	
	MeanCovarianceC GetZerothEstimate(const MeanCovarianceC &mc, const RealT &h, const RealT &dh, const UIntT &mval);
	//:Given a fin_pop estimate, this will obtain the T_0 and S_0 estimates from the data
	
	MeanCovarianceC FinMCEstimate(const MeanCovarianceC &mc, const MeanCovarianceC &mczero, const UIntT &m, const RealT &h);
	//:Compute least squares refined MCD estimate using Zeroth Step refinement
	
	RealT KSTest(const MeanCovarianceC &mc, const UIntT &m);
	//:Given a specific set of measurements, compute the K-S Test
private:
	UIntT imax;

};

#endif
