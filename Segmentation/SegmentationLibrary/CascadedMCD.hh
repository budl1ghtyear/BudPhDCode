#ifndef CascadedMCD_HH
#define CascadedMCD_HH
//////////////////////////////////////////////////////////////////////////////
//! file =      "Ravl/LipCode/Segmentation/SegmentationLibrary/CascadedMCD.hh"   
//! author =    "Bud Goswami"
//! lib =       SegLib
//! date =      "16.02.09"
//! userlevel = Normal
//! docentry =  "Ravl.API.LipCode.Segmentation.SegmentationLibrary"
//INCLUDE FILES
#include "MCD.hh"
using namespace RavlN;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------
// **********  CASCADEDMCD  **********************************************
// ---------------------------------------------------------------------------
//:Perform statistical estimation of data using greedy cluster extraction <p>
//:Class - declares the interface to the Cascaded Minimum Covariance Determinant Function<p>
//:Takes in as an input, a set of Sample of Vectors. <p>
//:Performs iterative MCD estimation on that set of vectors until it runs out of data<p>
//:Can return various final estimates from procedure.

class CascadedMCD : public MCD
{
	//Public Methods:
	//Default Constructor and Destructor
public:
	CascadedMCD() { cout<<"Empty CascadedMCD Object called"<<endl;}
	~CascadedMCD() {}
	
	
	CascadedMCD(const ImageC<RealRGBValueC> &rgbimg, const RealT &hval, const UIntT &numiter) : MCD(rgbimg,hval,numiter), minsamp(100){}
	//:Constructor inputs: ImageC<RealRGBValueC>, RealT (containing h) and UIntT (containing the number of iterations) <p>
	//variations of above (Constructor Must Always Have At Least ImageC<RealRGBValueC> object) <p>
	//More image variations will require recoding <p>
	//DEFAULTS VALUES:	<p>
	// h = 0.75 numiter = 10 minsamp =100 <p>
	CascadedMCD(const ImageC<RealRGBValueC> &rgbimg, const RealT &hval):MCD(rgbimg, hval), minsamp(100) {}
	CascadedMCD(const ImageC<RealRGBValueC> &rgbimg, const UIntT &numiter): MCD(rgbimg, numiter), minsamp(100){}
	CascadedMCD(const ImageC<RealRGBValueC> &rgbimg) : MCD(rgbimg), minsamp(100) {}
	
	CascadedMCD(const ImageC<VectorC> &vecimg, const RealT &hval, const UIntT &numiter) : MCD(vecimg, hval, numiter), minsamp(100) {}
	//:In case you want an image with just an arbitrary vector type
	CascadedMCD(const ImageC<VectorC> &vecimg, const RealT &hval) : MCD(vecimg,hval), minsamp(100) {}
	CascadedMCD(const ImageC<VectorC> &vecimg, const UIntT &numiter): MCD(vecimg, numiter), minsamp(100) {}
	CascadedMCD(const ImageC<VectorC> &vecimg) : MCD(vecimg), minsamp(100) {}	
	
	CascadedMCD(SampleData &smpl, const RealT &hval, const UIntT &numiter):MCD(smpl,hval), minsamp(100) {}
	//:Constructors with SampleC<VectorC>, RealT and UIntT <p>
	//:NOTE that using SampleDataC for images assumes that the vector includes index information <p>
	//:Form: vec[0] = imgdata0...vec[dim-1] = imgdata_dim; vec[dim] = Index.Row(), vec[dim+1] = Index.Col() <p>
	
	CascadedMCD(SampleData &smpl, const RealT &hval):MCD(smpl, hval), minsamp(100) {}
	CascadedMCD(SampleData &smpl, const UIntT &numiter) : MCD(smpl, numiter), minsamp(100) {}
	CascadedMCD(SampleData &smpl) : MCD(smpl), minsamp(100) {}
	
	CascadedMCD(const ImageC<TFVectorC<RealT,3> > &rgbimg, const RealT &hval, const UIntT &numiter):MCD(rgbimg,hval,numiter), minsamp(100) {}
	//:For usage with other 3D Real Colour Spaces
	CascadedMCD(const ImageC<TFVectorC<RealT,3> > &rgbimg, const RealT &hval):MCD(rgbimg,hval), minsamp(100) {}
	CascadedMCD(const ImageC<TFVectorC<RealT,3> > &rgbimg, const UIntT &numiter):MCD(rgbimg, numiter), minsamp(100) {}
	CascadedMCD(const ImageC<TFVectorC<RealT,3> > &rgbimg):MCD(rgbimg), minsamp(100) {}
	//Accessor Functions
	DListC< Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > Apply(void);
	UIntT GetMinSamples(void) const {return minsamp;}
	void SetMinSamples(UIntT &min) 
	{
		if(min >= 0) minsamp = min;
	 	else cerr<<"Min Samples must be a positive integer!";
	}
	void SetStartPop(const DListC<Tuple2C<VectorC,Index2dC> > &pop) 
	{
		start_pop.Empty();
		start_pop = pop.Copy();
	}
protected:
	//Private Data Members
	UIntT minsamp; //:to store minimum number of samples
};

#endif
