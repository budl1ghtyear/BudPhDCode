//PRecall.cc 
//Author: Bud Goswami
//Date: 16.01.08
//Description:
//Program calculates the precision recall measurements from a text file with segmentation quality scores

//Required Libs
#include "Ravl/Array1d.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/IO.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/ImageConv.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/IndexRange2d.hh"
#include "Ravl/Option.hh"
#include "Ravl/Stream.hh"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace RavlN;
using namespace RavlImageN;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//MAIN METHOD
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int nargs,char **argv) 
{
	//CREATE THE REQUIRED INPUT OPTIONS
	OptionC opt(nargs,argv);
	//IO Options
	FilenameC inf = opt.String("i","in.txt","Text File with Bounding Boxes");
	IntT numimg = opt.Int("n", 148, "Number of samples");
	RealT threshold = opt.Real("t",0.5,"Segmentation Quality Threshold");
	opt.Check();
	
	IStreamC is(inf);
	IntT truth = 0;
	//now we have the number of images, we run the loop numimg times
	for(IntT i = 0; i < numimg; i++)
	{
		RealT val = 0;
		is >> val;
		if(val >= threshold)
			truth ++;
	}
	RealT r = (RealT)truth/(RealT)numimg;
	cout<<"Ground Truth = "<<numimg<<"\n True = "<<truth<<"\t Threshold = "<<threshold<<"\t Recall = "<<r<<endl;
	return 0;
}
