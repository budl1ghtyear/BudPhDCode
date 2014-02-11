//LipClusteringC test program
//This is a new implementation to try and save a lot of time with the testing phase so that the algorithms only use the clustering step once
//Since the clustering step is going to be the most expensive, it should reduce the processing time immensely (by a factor of 3).

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
RealT ScoreImage(const IntT &gt, const ImageC<UIntT> &img, const FilenameC &img_file, const DirectoryC &g_dir);
ImageC<UIntT> MakeFullSizeImage(const ImageC<UIntT> &r_img, const ImageRectangleC &frame);
int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	DirectoryC img_dir = opt.String("i","Input/","example input image");
	IntT cluster_type = opt.Int("ct",0,"ClusteringType - values allowed - CMCD=0,KSMCD=1,KM=2,FCM=3, default = CMCD");
	IntT pixel_type = opt.Int("pt",0,"Pixel Type - RGB=0, NormRG=1, HSV=2, YUV=3, CIE=4, PHue=5 , default = RGB=0");
	//IntT label_type = opt.Int("lt",0,"Pixel labelling type - JNormal=0, JSkin=1, default =  JNormal=0 ");
	//IntT region_type = opt.Int("rt",0,"Region Identification type - CC=0,CCentral=1, default =CC=0 ");
	RealT hvalue = opt.Real("h",0.75,"Confidence Value of MCD - values between 0.5 and 1.0, default = 0.75");
	UIntT numclust = opt.Int("n",3,"Number of clusters for KMeans or Fuzzy C means, default = 3, must be greater than 1");
	UIntT gt_option = opt.Int("gt",0,"Ground Truth Type 0 - Basic Binary Image, 1 - AAM Data, 2 - Banca Data");
	DirectoryC gt_dir = opt.String("gd","GT/","Ground Truth Directory");
	opt.Compulsory("i");
	opt.Check();
	
	//Now instantiate the list of files that we need to start
	DListC<RealT> l0r0,l1r0,l0r1,l1r1; //to store the results of the computation
	DListC<StringC> img_files;
	switch(gt_option)
	{
		case 0:
		{
			img_files = img_dir.FiltList("*.ppm");
			break;
		}
		case 1:
		{
			img_files = img_dir.FiltList("*.ppm");
			break;
		}
		case 2:
		{
			img_files = img_dir.FiltList("*.png");
			break;
		}
	}
	//We have the file names of the files we are going to be processing, so just declare a loop to perform the processing
	//UIntT count = 0;
	for(DLIterC<StringC> it(img_files); it; it++)
	{
		//while(count < 10)
		{
		FilenameC img_name = img_dir + (*it);
		ImageC<RealRGBValueC> img;
		if(1==(!Load(img_name, img))) cerr<<"Could not load input image"<<endl;	
		//Just perform the modules one by one :) 
		LipSegmentationC lipseg(img_name);
		ColourSpaceImage col_img = lipseg.ColourSpaceConversion((PixelType)pixel_type);
		Tuple4C<Affine2dC, ImageRectangleC, ImageRectangleC,ImageC<VectorC> > amri = lipseg.AffineMouthRegionDetection(col_img);
		Tuple2C<MeanCovarianceC,MeanCovarianceC> col_trends = lipseg.LipClusteringStep((ClusteringType)cluster_type,ImageC<VectorC>(amri.Data4(),amri.Data2()),hvalue,numclust);
		//Now we should have the cluster data available to us, we can change the pixel labelling types and region identification types as we see fit
		//Perform scoring with l0 and r0
		ImageC<UIntT> ppl_img = lipseg.ProbabilisticPixelLabellingStep(ImageC<VectorC>(col_img.GetImage(),amri.Data3()) , col_trends.Data1(), col_trends.Data2(), (LabellingType)0);
		ImageC<UIntT> rgn_img = lipseg.RegionIdentificationStep(ppl_img,(RegionIdentificationType)0);
		ImageC<UIntT> lip_img = MakeFullSizeImage(rgn_img, img.Frame());
		//if(!Save("@X:L0R0 Image",lip_img)) cerr<<"Could not show saved image"<<endl;
		l0r0.Append(ScoreImage(gt_option, lip_img, img_name, gt_dir)); 
		//Perform scoring with l1 and r0
		ppl_img = lipseg.ProbabilisticPixelLabellingStep(ImageC<VectorC>(col_img.GetImage(),amri.Data3()) , col_trends.Data1(), col_trends.Data2(), (LabellingType)1);
		rgn_img = lipseg.RegionIdentificationStep(ppl_img,(RegionIdentificationType)0);
		lip_img = MakeFullSizeImage(rgn_img, img.Frame());
		//if(!Save("@X:L1R0 Image",lip_img)) cerr<<"Could not show saved image"<<endl;
		l1r0.Append(ScoreImage(gt_option, lip_img, img_name, gt_dir)); 
		//Perform scoring with l0 and r1
		ppl_img = lipseg.ProbabilisticPixelLabellingStep(ImageC<VectorC>(col_img.GetImage(),amri.Data3()) , col_trends.Data1(), col_trends.Data2(), (LabellingType)0);
		rgn_img = lipseg.RegionIdentificationStep(ppl_img,(RegionIdentificationType)1);
		lip_img = MakeFullSizeImage(rgn_img, img.Frame());
		//if(!Save("@X:L0R1 Image",lip_img)) cerr<<"Could not show saved image"<<endl;
		l0r1.Append(ScoreImage(gt_option, lip_img, img_name, gt_dir)); 
		//Perform scoring with l1 and r1
		ppl_img = lipseg.ProbabilisticPixelLabellingStep(ImageC<VectorC>(col_img.GetImage(),amri.Data3()) , col_trends.Data1(), col_trends.Data2(), (LabellingType)1);
		rgn_img = lipseg.RegionIdentificationStep(ppl_img,(RegionIdentificationType)1);
		lip_img = MakeFullSizeImage(rgn_img, img.Frame());
		//if(!Save("@X:L1R1 Image",lip_img)) cerr<<"Could not show saved image"<<endl;
		l1r1.Append(ScoreImage(gt_option, lip_img, img_name, gt_dir)); 
		//count++;
		}
	}
	//Now we should have a complete set of DLists with the various lip pixel labelling and region identification options, just print them out
	cout<<" L0R0 \n" <<l0r0<<endl;
	cout<<"\n L1R0 \n" <<l1r0<<endl;
	cout<<"\n L0R1 \n" <<l0r1<<endl;
	cout<<"\n L1R1 \n" <<l1r1<<endl;
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

RealT ScoreImage(const IntT &gt, const ImageC<UIntT> &img, const FilenameC &img_file, const DirectoryC &g_dir)
{
	RealT score = 0.0;
	//Grab the ground truth
	switch(gt)
	{
		//Case for regular manually annotated ground truth
		case 0:
		{
			FilenameC gt_file = g_dir + img_file.BaseNameComponent() + ".pbm";
			ImageC<UIntT> gt_img;
			if(!Load(gt_file,gt_img)) cerr<<"Ground Truth Image Could Not Be Loaded"<<endl;
			score = ComputeSegmentationQuality(img,gt_img);
			break;
		}
		case 1:
		{	
			//Jean Yves AAM Markup Data
			FilenameC gt_file = g_dir + img_file.BaseNameComponent() + ".xml";
			XMLToGTConv gt_conv(gt_file, img.Frame());
			ImageC<UIntT> gt_img = gt_conv.Apply().Copy();
			//cout << gt_img.Frame() << " " << lip_img.Frame() << endl;
			score = ComputeSegmentationQuality(img,gt_img);
			break;
		}
		case 2:
		{
			//BANCA Data
			FilenameC gt_file = g_dir + img_file.BaseNameComponent() + ".did";
			//cout<<"GT Filename is - "<<gt_file<<endl;
			BANCAGT gt_conv(gt_file,img.Frame());
			ImageC<UIntT> gt_img = gt_conv.GetGTImage();
			//if(!Save("@X:BANCA GT IMage",gt_img)) cerr<<"Could not save BANCA GT Image"<<endl;
			score = ComputeSegmentationQuality(img,gt_img);
			break;
		}
		default:
		{
			break;
		}
	}
	return score;		
}

ImageC<UIntT> MakeFullSizeImage(const ImageC<UIntT> &r_img, const ImageRectangleC &frame)
{
	ImageC<UIntT> res(frame,0);
	for(Array2dIterC<UIntT> it(r_img); it; it++)
	{
		if((*it) == 255)
			res[it.Index()] = 255;
	}	
	return res;
} 
