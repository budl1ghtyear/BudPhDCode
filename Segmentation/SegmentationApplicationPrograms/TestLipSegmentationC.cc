//LipClusteringC test program

//Just to test the Lip Clustering Test Class
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/Option.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/String.hh"
#include "LipSegmentation.hh"
#include "ColourSpaceImage.hh"
#include "XMLToGTConv.hh"
#include "BANCAGT.hh"
#include "CostFunctions.hh"
using namespace RavlN;
using namespace RavlImageN;

//RealT ComputeSegmentationQuality(const ImageC<UIntT> &im, const ImageC<UIntT> &gt);
//Method to compute the image quality based on binary ground truth and segmented lip image

int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	FilenameC img_name = opt.String("i","Input.ppm","example input image");
	IntT cluster_type = opt.Int("ct",0,"ClusteringType - values allowed - CMCD=0,KSMCD=1,KM=2,FCM=3, default = CMCD");
	IntT pixel_type = opt.Int("pt",0,"Pixel Type - RGB=0, NormRG=1, HSV=2, YUV=3, CIE=4, PHue=5 , default = RGB=0");
	IntT label_type = opt.Int("lt",0,"Pixel labelling type - JNormal=0, JSkin=1, default =  JNormal=0 ");
	IntT region_type = opt.Int("rt",0,"Region Identification type - CC=0,CCentral=1, default =CC=0 ");
	RealT hvalue = opt.Real("h",0.75,"Confidence Value of MCD - values between 0.5 and 1.0, default = 0.75");
	UIntT numclust = opt.Int("n",3,"Number of clusters for KMeans or Fuzzy C means, default = 3, must be greater than 1");
	UIntT gt_option = opt.Int("gt",0,"Ground Truth Type 0 - Basic Binary Image, 1 - AAM Data, 2 - Banca Data");
	DirectoryC gt_dir = opt.String("gd","GT/","Ground Truth Directory");
	opt.Compulsory("i");
	opt.Check();
	
	//Load input image
	ImageC<RealRGBValueC> img;
	if(1==(!Load(img_name, img))) cerr<<"Could not load input image"<<endl;
	//if(!Save("@X:Original Image",img)) cerr<<"Could not show input image"<<endl;
	//cout<<"Instantiating class"<<endl;
	LipSegmentationC lipseg(img_name);
	Tuple3C<ImageC<UIntT>, MeanCovarianceC, MeanCovarianceC> result = lipseg.Apply((PixelType)pixel_type,(ClusteringType)cluster_type,(LabellingType) label_type,(RegionIdentificationType) region_type,hvalue, numclust);
	//if(!Save("@X:Binary Image",result.Data1())) cerr<<"Could not show binary lip file"<<endl; //Now we have the lip binary image in real space
	//cout<<"Lip FOund"<<endl;
	//Make full size Binary Lip Image
	ImageC<UIntT> lip_img(img.Frame(),0);
	for(Array2dIterC<UIntT> it(result.Data1()); it; it++)
	{
		if((*it) == 255)
			lip_img[it.Index()] = 255;
	}
	RealT score = 0.0;
	//Grab the ground truth
	switch(gt_option)
	{
		//Case for regular manually annotated ground truth
		case 0:
		{
			FilenameC gt_file = gt_dir + img_name.BaseNameComponent() + ".pbm";
			ImageC<UIntT> gt_img;
			if(!Load(gt_file,gt_img)) cerr<<"Ground Truth Image Could Not Be Loaded"<<endl;
			score = ComputeSegmentationQuality(lip_img,gt_img);
			break;
		}
		case 1:
		{	
			//Jean Yves AAM Markup Data
			FilenameC gt_file = gt_dir + img_name.BaseNameComponent() + ".xml";
			XMLToGTConv gt_conv(gt_file, lip_img.Frame());
			ImageC<UIntT> gt_img = gt_conv.Apply().Copy();
			//cout << gt_img.Frame() << " " << lip_img.Frame() << endl;
			score = ComputeSegmentationQuality(lip_img,gt_img);
			break;
		}
		case 2:
		{
			//BANCA Data
			FilenameC gt_file = gt_dir + img_name.BaseNameComponent() + ".did";
			//cout<<"GT Filename is - "<<gt_file<<endl;
			BANCAGT gt_conv(gt_file,lip_img.Frame());
			ImageC<UIntT> gt_img = gt_conv.GetGTImage();
			score = ComputeSegmentationQuality(lip_img,gt_img);
			break;
		}
		default:
		{
			break;
		}
	}
	cout<<score<<endl;
	return 0;
}
//~ 
//~ RealT ComputeSegmentationQuality(const ImageC<UIntT> &im, const ImageC<UIntT> &gt)
//~ {
	//~ RealT num = 0.0,den=0.0;
	//~ for(Array2dIter2C<UIntT, UIntT> it(im,gt); it; it++)
	//~ {
		//~ if(it.Data1() == 255)
		//~ {
			//~ den++;
			//~ if(it.Data2() == 255)
			//~ {
				//~ num++;
			//~ }
		//~ }
		//~ else
		//~ {
			 //~ if(it.Data2() == 255)
			//~ {
				//~ den++;
			//~ }
		//~ }
	//~ }
	//~ return ((RealT)(num/den));
//~ }
