//This is the zeroth order PCAStateSpace tracker that uses the active contour measurement model

//INCLUDE FILES:
#include "Ravl/Option.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/DP/SequenceIO.hh"
#include "Ravl/Image/DeinterlaceStream.hh"

#include "BasisSpline.hh"
#include "Particle.hh"
#include "SystematicResampling.hh"
#include "SIR.hh"
#include "LipSegmentation.hh"
#include "ColourSpaceImage.hh"
#include "XMLToGTConv.hh"
#include "BANCAGT.hh"
#include "LipExtractorC.hh"
#include "TrackingUtilityFunctions.hh"
#include "PCAStateVector.hh"
#include "ZOAVPropModel.hh"
#include "ACMeasurementModel.hh"
#include "EdgeMeasurementModel.hh"
#include "ICAStateVector.hh"
#include "BANCAGT.hh"
#include "Ravl/Assert.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Image/DeinterlaceStream.hh"
using namespace RavlN;
using namespace RavlImageN;
RandomGaussC BasePropagationModelC::rnd;

void DrawCrosses(ImageC<RealRGBValueC> &im, SArray1dC<Point2dC> &pts);

int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	FilenameC vid_file = opt.String("i","./videofile.avi","Input Video File");
	FilenameC first_file = opt.String("f","./firstfile.did","The coordinates of the first lip shape that was bootstrapped. This is in .did format");
	DirectoryC pca_dir = opt.String("p","/vol/vssp/lip-tracking/StatisticalShapeModelsOfLips/PCA/","PCA Data Directory");
	UIntT subspace_proj = opt.Int("s",0,"Method of subspace projection to use- 0: PCA, 1:ICA");
	UIntT numcomponents = opt.Int("npca",10,"Number of eigen-components to choose in state representation");
	//Particle Filter Parameters
	//~ UIntT search = opt.Int("s",20,"Length of search window in pixels in affine space");
	//~ UIntT noise = opt.Int("n",2,"Length of noise perturbtaion in pixels");
	//~ IntT numpts = opt.Int("pt",50,"Number of Particles");
	opt.Compulsory("i");
	opt.Check();
	
	DPIPortC<ImageC<RealRGBValueC> > in;
	RavlAlwaysAssertMsg(!OpenISequence(in, vid_file),StringC("Could not open ") + vid_file + " for i/p"); 
	cout<<"Loaded video stream"<<endl;
	DeinterlaceStreamC<RealRGBValueC> din(in);
	ImageC<RealRGBValueC> im;
	cout<<"Loaded deinterlaced stream"<<endl;
	UIntT i=1;
	SArray1dC<Point2dC> first_pts;
	cout<<"Start frame: "<<din.Start()<<endl;
	while(din.GetAt(i,im))
	{
		if(!Save("@X:Img File", im)) exit(1);
		if(i == 1)
		{
			//~ BANCAGT bancabootstrap(first_file,im.Frame());
			//~ DListC<Point2dC> banca_outerpts = bancabootstrap.GetOuterPoints();
			//~ SArray1dC<Point2dC> dummy(banca_outerpts.Size());
			//~ for(UIntT pt_ind = 0; pt_ind < banca_outerpts.Size(); pt_ind++)
			//~ {
				//~ dummy[pt_ind] = banca_outerpts.Nth(pt_ind).Copy();
			//~ }
			//~ first_pts = dummy.Copy();
		}
		
		cout<<"frame idx: "<<i<<" "<<din.Tell64()<<endl;
		i+=2;
	}
	cout<<"Size: "<<din.Size()/2<<endl;
	DrawCrosses(im, first_pts);
	if(!Save("@X:First Lips with Img File", im)) exit(1);
	return 0;
}

void DrawCrosses(ImageC<RealRGBValueC> &im, SArray1dC<Point2dC> &pts)
{
	RealRGBValueC col(255,10,10);
	for(SArray1dIterC<Point2dC> it(pts); it; it++)
	{
		Index2dC ind((*it).Row(),(*it).Col());
		DrawCross(im,col,ind,4);
	}
	if(!Save("@X:Cross Image",im)) cerr<<"Could not display crosses image"<<endl;
}
