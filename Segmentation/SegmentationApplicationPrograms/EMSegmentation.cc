//File: EMSegmentation.cc
//Input: Query Image Filename, Ground Truth Image Path, Number of Clusters
//Output: Displays the segmented lip region
//Author: Bud Goswami
//Date: 17.02.09
//Data Structures and Iterators
#include "Ravl/Array1d.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/Array2dIter2.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/Tuple2.hh"
//OS Stuff
#include "Ravl/IO.hh"
#include "Ravl/Option.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/OS/Filename.hh"
//Image Stuff
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/ImageConv.hh"
#include "Ravl/Image/DrawFrame.hh"
#include "Ravl/Image/OpenCVConvert.hh"
#include "Ravl/Image/ImageRectangle.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Ellipse2d.hh"
//Pattern Recognition Stuff
#include "Ravl/MeanCovariance.hh"
#include "Ravl/Image/ConnectedComponents.hh"
#include "Ravl/Image/Segmentation.hh"
#include "Ravl/PatternRec/DesignGaussianMixture.hh"
#include "Ravl/PatternRec/Sample.hh"
#include "Ravl/Vector.hh"
#include "Ravl/SArray1d.hh"
#include "Ravl/PatternRec/Function.hh"
#include "Ravl/PatternRec/DesignKMeans.hh"
#include "Ravl/PatternRec/DistanceMahalanobis.hh"
#include "Ravl/PatternRec/Classifier.hh"
#include "Ravl/PatternRec/DesignClassifierGaussianMixture.hh"
using namespace RavlImageN;

int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	FilenameC qimg = opt.String("i","in.ppm","Input Query Image File");
	//DirectoryC gtdir = opt.String("g","gt.pbm","Directory Containing Ground Truth Image");
	//MCD Related Parameters
	UIntT numclust = opt.Int("c",3,"Expected Number of Clusters");
	opt.Check();
	
	//Load Image
	ImageC<RealRGBValueC> src;   
	if(!Load(qimg, src)) cerr<<"Loading RAVL Image Failed"<<endl;
	
	//Insert Image Data into SampleC<VectorC> after inserting into SArray1dC
	//Using RGB Colour Representation
	SArray1dC<VectorC> dummy(src.Frame().Area());
	UIntT ind = 0;
	for(Array2dIterC<RealRGBValueC> it(src); it; it++)
	{
		VectorC v((*it).Red(),(*it).Green(),(*it).Blue());
		dummy[ind] = v.Copy();
		ind++;
	}
	SArray1dC<UIntT> outdummy(src.Frame().Area()); outdummy.Fill(0);
	cout<<"Calculated Sample VectorC"<<endl; //dummy<<endl;
	SampleC<VectorC> sample(dummy);
	SampleC<UIntT> output(outdummy);
	UIntT mixtures = 3;
	//Declare new GaussianMixture Classifier
	DesignClassifierGaussianMixtureC emcl(3);
	//Create classifier class
	ClassifierC cl = emcl.Apply(sample, output);
	//cout<<output<<endl;
	ImageC<UIntT> clustimg(src.Frame(),0);
	//Perform classification given designed classifier on the image data.
	for(Array2dIterC<RealRGBValueC> it(src); it; it++)
	{
		VectorC v((*it).Red(), (*it).Green(),(*it).Blue());
		UIntT label = cl.Classify(v);
		clustimg[it.Index()] = label;
	}
	SegmentationC seg(clustimg, 3);
 	if(!Save("@X: Clustered Image",seg.RandomImage())) cerr<<"Could not show file"<<endl;
	return 1;
	
}
