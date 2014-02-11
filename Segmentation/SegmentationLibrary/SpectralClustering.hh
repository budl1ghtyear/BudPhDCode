//File: SpectralClustering.cc

//Date: 17.02.09

//Datavariables:
//Private Data Members
//Algorithm Parameters


#ifndef SpectralClustering_HH
#define SpectralClustering_HH
//////////////////////////////////////////////////////////////////////////////
//! file =      "Ravl/LipCode/Segmentation/SegmentationLibrary/SpectralClustering.hh"   
//! author =    "Bud Goswami"
//! lib =       SegLib
//! date =      "17.02.09"
//! userlevel = Normal
//! docentry =  "Ravl.API.LipCode.Segmentation.SegmentationLibrary"
//INCLUDE FILES

//Data Structures & Iterators
#include "Ravl/DArray1dIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Matrix.hh"
#include "Ravl/Matrix2d.hh"
#include "Ravl/Image/Image.hh"
#include "SampleData.hh"
#include "Ravl/Tuple2.hh"
#include "Ravl/Tuple3.hh"
#include "Ravl/Vector.hh"
//Pixel Types
#include "Ravl/Image/RealRGBValue.hh"
//Pattern Recognition Classes
#include "Ravl/PatternRec/Classifier.hh"
#include "Ravl/PatternRec/DesignFuncPCA.hh"
#include "Ravl/PatternRec/DesignKMeans.hh"
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
// **********  Spectral Clustering  **********************************************
// ---------------------------------------------------------------------------
//:Class - defines the interface to the Spectral Clustering Function <p>
//<p>Takes in as an input, a set of Sample of Vectors. Performs MCD estimation on that set of vectors.</p>
//Can return various final estimates from procedure. <p>
// For more details refer to U Von Luxberg (excellent video tutorial on videolectures.net)
class SpectralClustering
{
	enum SimType{GAUSSIANKERNEL, SHINCUT};
	//: Enumeration for the similarity matrix generation method
	enum LapType{RANDOMWALK, UNNORMALISED, SYMMETRIC};
	//: Enumeration for the Laplacian matrix generation method

	//Public Methods:
	//Default Constructor and Destructor
public:
	SpectralClustering() { cout<<"Empty MCD Object called"<<endl;}
	~SpectralClustering() {}
	
	
	SpectralClustering(const ImageC<RealRGBValueC> &rgbimg, const UIntT &cl);
	//:Constructor inputs: ImageC<RealRGBValueC>, RealT (containing h) and UIntT (containing the number of iterations)
	SpectralClustering(const ImageC<RealRGBValueC> &rgbimg);
	SpectralClustering(const ImageC<VectorC> &vecimg, const UIntT &cl);
	SpectralClustering(const ImageC<VectorC> &vecimg);

	//Accessor Functions
	SampleC<VectorC> GetObservationSamples(void) const {return obs;}
	void SetObservationSamples(SampleC<VectorC> &sm)
	{
		obs = sm.Copy();
	}
	UIntT GetNumClusters(void) const {return numclust;}	
	void SetNumClusters(const UIntT &nc) 
	{
		if(nc > 1)
			numclust = nc;
		else
		{
			cerr<<"Number of clusters MUST be greater than 1, incorporating 1 eigvect for 0 answer"<<endl;
			numclust = 3;
		}
	}

	//Public Method
	Tuple2C<MatrixC,MatrixC> ComputeSDMatrix(SimType s);
	//:Compute the DS Matrix (notation from Luxburg paper)
	MatrixC ComputeLaplacian(const Tuple2C<MatrixC, MatrixC> &sd, LapType L);
	//: Compute the laplacian matrix
	MatrixC GetNormalisedEigenVectors(MatrixC &lmat);
	//: GetNormalised Eigen Vectors
	DListC<UIntT> GetKMeansClassification(const MatrixC &emat);
	//: Get the KMeans Classification results
	
	ImageC<UIntT> Apply();
	//:Apply method
	//:Apply method
protected:
	//Data members
	SampleC<VectorC> obs; //:to store the initial observations
	DListC<Tuple2C<VectorC, Index2dC> > obsind; //:to store relationship between vector and location
	UIntT numclust; //:store the number of clusters
	ImageRectangleC imrect;//: To store the ROI so we only perform the operation on a small amount of data.
	//Functions
	Tuple2C<SampleC<VectorC>, DListC<Tuple2C<VectorC, Index2dC> > > ConvertImageToSamples(const ImageC<RealRGBValueC> &rgb);
	//: Data structure conversion functions
	Tuple2C<SampleC<VectorC>, DListC<Tuple2C<VectorC, Index2dC> > > ConvertImageToSamples(const ImageC<VectorC> &rgb);
	//: Data structure conversion functions
	
};

#endif
