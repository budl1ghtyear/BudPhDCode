//File: MCD.hh
//Class - declares the interface to the Minimum Covariance Determinant Function
// Takes in as an input, a set of Sample of Vectors. Performs MCD estimation on that set of vectors.
// Can return various final estimates from procedure. 
//Author: Bud Goswami
//Date: 29.01.09

#ifndef MCD_HH
#define MCD_HH
//////////////////////////////////////////////////////////////////////////////
//! file =      "Ravl/LipCode/Segmentation/SegmentationLibrary/MCD.hh"   
//! author =    "Bud Goswami"
//! lib =       Bspline
//! date =      "29/01/09"
//! userlevel = Normal
//! docentry =  "Ravl.API.LipCode.Segmentation.SegmentationLibrary"
//INCLUDE FILES

//Data Structures & Iterators
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Image/Image.hh"
#include "SampleData.hh"
#include "Ravl/TFVector.hh"
#include "Ravl/Tuple2.hh"
#include "Ravl/Tuple3.hh"
#include "Ravl/Vector.hh"
//Pixel Types
#include "Ravl/Image/RealRGBValue.hh"
//Pattern Recognition Classes
#include "Ravl/MeanCovariance.hh"
#include "Ravl/PatternRec/Sample.hh"
#include "Ravl/Random.hh"
#include "Ravl/SumsNd2.hh"
//Misc Classes
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/ImageRectangle.hh"

using namespace RavlN;
using namespace RavlImageN;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------
// **********  MCD  **********************************************
// ---------------------------------------------------------------------------
// (previous lines will be ignored - no preceeding //: )

//: Class Implementing Minimum Covariance Determinant Method for Robust Statistical Estimation <p>
//
// This algorithm implements FASTMCD for robust statistical estimation of data as specified in the paper by Rousseuw

class MCD
{
	//Public Methods:
	 
public:
	MCD() { cout<<"Empty MCD Object called"<<endl;}
	//:Default Constructor
	~MCD() {}
	//:Destructor
	
	MCD(const ImageC<RealRGBValueC> &rgbimg, const RealT &hval, const UIntT &numiter);
	//:Constructor inputs: ImageC<RealRGBValueC>, RealT (containing h) and UIntT (containing the number of iterations) <p>
	//:variations of above (Constructor Must Always Have At Least ImageC<RealRGBValueC> object)<p>
	//:More image variations will require recoding<p>
	//:DEFAULTS VALUES:	<p>
	//:h = 0.75 numiter = 10;
	MCD(const ImageC<RealRGBValueC> &rgbimg, const RealT &hval);
	MCD(const ImageC<RealRGBValueC> &rgbimg, const UIntT &numiter);
	MCD(const ImageC<RealRGBValueC> &rgbimg);
	
	MCD(const ImageC<VectorC> &vecimg, const RealT &hval, const UIntT &numiter);
	//:In case you want an image with just an arbitrary vector type
	MCD(const ImageC<VectorC> &vecimg, const RealT &hval);
	MCD(const ImageC<VectorC> &vecimg, const UIntT &numiter);
	MCD(const ImageC<VectorC> &vecimg);	
	
	MCD(SampleData &smpl, const RealT &hval, const UIntT &numiter);
	//:Constructors with SampleC<VectorC>, RealT and UIntT <p>
	//:NOTE that using SampleDataC for images assumes that the vector includes index information <p>
	//:Form: vec[0] = imgdata0...vec[dim-1] = imgdata_dim; vec[dim] = Index.Row(), vec[dim+1] = Index.Col() <p>
		
	MCD(SampleData &smpl, const RealT &hval);
	MCD(SampleData &smpl, const UIntT &numiter);
	MCD(SampleData &smpl);
	
	MCD(const ImageC<TFVectorC<RealT,3> > &rgbimg, const RealT &hval, const UIntT &numiter);
	//:Additional constructors that take in any other RAVL 3D pixel types
	MCD(const ImageC<TFVectorC<RealT,3> > &rgbimg, const RealT &hval);
	MCD(const ImageC<TFVectorC<RealT,3> > &rgbimg, const UIntT &numiter);
	MCD(const ImageC<TFVectorC<RealT,3> > &rgbimg);
	
	
	//:Accessor Functions
	RealT GetH(void) {return hvalue;}
	void SetH(const RealT &h)
	{
		if((h > 0.5)&&(h <= 1.0))
			hvalue = h;
		else
		{ 
			cerr<<"H Range is 0.5<h<1.0, setting to default - 0.75"<<endl;
			hvalue = 0.75;
		}
	}
	UIntT GetN(void) {return nvalue;}
	DListC<Tuple2C<VectorC,Index2dC> > GetStartPop(void) {return start_pop;}
	Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > GetElementalSet(void);
	Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > Apply(void);
	
protected:
	//Private Data Members
	//Algorithm Parameters
	RealT hvalue; //:to store confidence in parameters
	UIntT nvalue; //:to store number of iterations to convergence
	DListC<Tuple2C<VectorC,Index2dC> > start_pop;
	//:Starting data population: Form(DListC<Tuple3C<VectorC,Index2dC>>) so each data vector has an allocated index
	
	Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > elemental_set;
	//:Elemental Set: Form Tuple2C<MeanCovarianceC,DListC<Tuple3C<VectorC,Index2dC> > >
	
	Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > end_pop;
	//:Final Set: Same as above
	//Private Member Functions
	//:Data Converter Functions to ensure input args are conformant to start_pop expectations
	DListC<Tuple2C<VectorC,Index2dC> > GetStartPopFromImage(const ImageC<RealRGBValueC> &img); 
	DListC<Tuple2C<VectorC,Index2dC> > GetStartPopFromSampleData(SampleData &img);
	DListC<Tuple2C<VectorC,Index2dC> > GetStartPopFromImageVector(const ImageC<VectorC> &img);
	DListC<Tuple2C<VectorC,Index2dC> > GetStartPopFromImageVector(const ImageC<TFVectorC<RealT,3> > &img);
	
	MeanCovarianceC GetCStep(const MeanCovarianceC &mcstart, const RealT &h, const  DListC<Tuple2C<VectorC,Index2dC> > &st_pop);
	//:Perform C-Step
};

#endif
