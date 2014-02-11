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
using namespace RavlN;
using namespace RavlImageN;

int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	DirectoryC img_dir = opt.String("i","Input/","example input image");
	IntT cluster_type = opt.Int("ct",0,"ClusteringType - values allowed - CMCD=0,KSMCD=1,KM=2,FCM=3, default = CMCD");
	IntT pixel_type = opt.Int("pt",0,"Pixel Type - RGB=0, NormRG=1, HSV=2, YUV=3, CIE=4, PHue=5 , default = RGB=0");
	UIntT gt_option = opt.Int("gt",0,"Ground Truth Type 0 - Basic Binary Image, 1 - AAM Data, 2 - Banca Data");
	DirectoryC gt_dir = opt.String("gd","GT/","Ground Truth Directory");
	opt.Compulsory("i");
	opt.Check();
	
	//Now instantiate the list of files that we need to start
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
	for(DLIterC<StringC> it(img_files); it; it++)
	{
		FilenameC img_name = img_dir + (*it);
		ImageC<RealRGBValueC> img;
		if(1==(!Load(img_name, img))) cerr<<"Could not load input image"<<endl;	
		//Compute the Mouth Region Image
		FilenameC  faceDetectionModelHQ = "/vol/vssp/localsoft/External/FaceDetect/models/HighQualityFD.abs";
		FilenameC  faceDetectionModelLQ = "/vol/vssp/localsoft/External/FaceDetect/models/LowQualityFaceScan.abs";	
		FaceScanC detectAAMFace;           // face locator that can provide face shape from AAM
		DetectFaceC detectFace;            // eezy-use face locator
		DetectedFaceC detectedFace;
		//cout<<"set up face and eye locators"<<endl;
		FaceDetectionModelC faceDetectionModel;
		if(!Load(faceDetectionModelHQ, faceDetectionModel))
			RavlIssueError("Trouble loading eezy-use face detection model\n");
		detectFace = DetectFaceC(faceDetectionModel); 
		if(!Load(faceDetectionModelLQ, detectAAMFace))
			RavlIssueError("Trouble loading AAM face detection model\n");
		ImageC<ByteRGBValueC> fullSizeImg = RealRGBImageCT2ByteRGBImageCT(img);
		UIntT factor(1);
		//cout<<"DetectFace.Apply()"<<endl;
		DetectedFaceC tmp = detectFace.Apply(fullSizeImg);
		detectedFace = DetectedFaceC(fullSizeImg, tmp.LeftEye(), tmp.RightEye(), tmp.Quality());
		//cout<<"Convert to output format"<<endl;
 
}	

	}
	return 0;
}


