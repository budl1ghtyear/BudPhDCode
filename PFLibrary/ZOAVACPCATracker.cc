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
using namespace RavlN;
using namespace RavlImageN;
RandomGaussC BasePropagationModelC::rnd;

void DrawCrosses(ImageC<RealRGBValueC> &im, SArray1dC<Point2dC> &pts);

int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	DirectoryC img_dir = opt.String("i","./Images/","Input Image Directory");
	DirectoryC pca_dir = opt.String("p","/vol/vssp/lip-tracking/StatisticalShapeModelsOfLips/PCA/","PCA Data Directory");
	//DirectoryC pca_dir = opt.String("p","/vol/vssp/lip-tracking/StatisticalShapeModelsOfLips/JADE/RedJADE/","PCA Data Directory");
	StringC file_ext = opt.String("e","ppm","QUERY INPUT FILE EXTENSION");
	//PCA Parameters
	UIntT numcomponents = opt.Int("npca",10,"Number of eigen-components to choose in state representation");
	//Particle Filter Parameters
	UIntT search = opt.Int("s",20,"Length of search window in pixels in affine space");
	UIntT noise = opt.Int("n",2,"Length of noise perturbtaion in pixels");
	IntT numpts = opt.Int("pt",50,"Number of Particles");
	opt.Compulsory("i");
	opt.Check();

	DListC<StringC> img_files = img_dir.FiltList("*."+file_ext); 
	FilenameC first_file = img_dir + img_files.First();
	DPIPortC<ImageC<RealRGBValueC> > inputstream;
	if(!OpenISequence(inputstream, (StringC)first_file))
	{
		cerr<<"Failed to open input sequence."<< img_dir<<" \n " ;
		return 1;
	}
	DeinterlaceStreamC<RealRGBValueC> dinputstream(inputstream);
	ImageC<RealRGBValueC> first_img = dinputstream.Get();
	RealRGBValueC color(10,200,10); RealRGBValueC color2(10,10,200);
	cout<<"Instantiated Streams"<<endl;
	///////////////////////////
	//BOOTSTRAP PROCESS
	//Load input image
	Tuple2C<Point2dC,Point2dC> first_eyes;
	//Loop until we get a successful eye detection, the process starts at this point
	while(!dinputstream.IsGetEOS())
	{
		try
		{
			first_eyes= GetEyeLoc(first_img);
		}
		catch(OmniN::ColossusExceptionNoFaceFoundC er)
		{
			//MOVE ON TO THE NEXT FRAME		
			first_img = dinputstream.Get();
			continue;
		}
		break;
	}
	//At this point, we have successfully detected a set of eyes
	Tuple2C<BasisSplineC,Affine2dC> bootstrap = GetLipPoints(first_img);
	//Now affine transform the lip points onto face space
	SArray1dC<Point2dC> aff_spl = AffineTransformPoints(bootstrap.Data1().GetControlPoints(),bootstrap.Data2());
	Tuple2C<ImageC<RealRGBValueC>,Affine2dC> aff_tr = AffNormWRTEyes(first_img,first_eyes.Data1(),first_eyes.Data2());
	//Draw the resulting points on top of the affine transformed image in "head-space";
	ImageC<RealRGBValueC> headim(aff_tr.Data1().Copy());
	RealRGBValueC cross(255,10,10);
	if(!Save("@X:Displaying affine transformed head image",headim)) cerr<<"Could not show the head image"<<endl;
	for(SArray1dIterC<Point2dC> it(aff_spl); it; it++)
		DrawCross(headim,cross,(*it));
	if(!Save("@X:Displaying affine transformed head points",headim)) cerr<<"Could not show the head image"<<endl;
	cout<<"Aff_Spl = "<<aff_spl<<endl;
	BasisSplineC affspl(aff_spl);
	//Now we are ready to instantiate the Particle Filter:
	//PCAStateVectorC initpts(pca_dir, affspl.GetControlPoints(),numcomponents);
	ICAStateVectorC initpts(pca_dir, affspl.GetControlPoints());
	RealT initwt = 1.0 / (RealT) numpts;
	ParticleC initparticle(initpts,initwt);
	Array1dC<ParticleC> ptcol(numpts);
	ptcol.Fill(initparticle);	
	cout<<"Initial Particle state "<<initparticle.GetState().GetStateVector()<<endl;
	cout<<"Initial Particle weight "<<initparticle.GetWeight()<<endl;
	UIntT nth = 1;
	SArray1dC<Point2dC> prev_snaxels = aff_spl.Copy(); //This contains the lip B-Spline control points in affine transformed space
	//TRACKING LOOP INSTANTIATES HERE:
	while(!dinputstream.IsGetEOS())
	{
		ImageC<RealRGBValueC> img = dinputstream.Get();
		Tuple2C<Point2dC,Point2dC> img_eyes;
		try
		{
			img_eyes = GetEyeLoc(img);			
		}
		// IF THE EYE DETECTION FAILS, THEN WE INSTANTIATE IT TO THE PREVIOUS EYE LOCATION 
		// THIS WILL FAIL IF THE HEAD MOTION FOR THIS FRAME IS HUGELY DIFFERENT TO THE PREVIOUS FRAME! 
		catch(OmniN::ColossusExceptionNoFaceFoundC er)
		{
			img_eyes = first_eyes.Copy();
			continue;
		}
		//NOW AFFINE TRANSFORM THE TRACKING IMAGE
		Tuple2C<ImageC<RealRGBValueC>,Affine2dC> aff_tr = AffNormWRTEyes(img,img_eyes.Data1(),img_eyes.Data2());
		//Instantiate the propagation model
		ZOAVPropModelC ofprop(noise);
		//Instantiate the measurement model
		//~ Tuple2C<ImageC<RealT>, ImageC<Vector2dC> > measims = GetMeasurementImages(img,img_eyes,first_img,first_eyes);
		//~ if(!Save("@X:Susan Edge Image ",measims.Data1())) cerr<<"Could not load the edge image produced by Susan"<<endl;
		//~ ACMeasurementModelC meas(measims.Data1(),measims.Data2(),prev_snaxels);
		EdgeMeasurementModelC meas(aff_tr.Data1());
		//Instantiate the resampling model
		SystematicResamplingC sysres;
		//Instantiate the particle filter
		SIR pfilter(numpts, ofprop, meas, sysres,ptcol);
		cout<<"Starting Particle Filter here"<<endl;
		Array1dC<Point2dC> pf_pts = pfilter.Apply().GetObservationVector(); //GET THE APPLY METHOD 
		//DisplayImage(ArrayToSArray_B(pf_pts));
		//Now we have a result for the Particle Filter Estimate, just obtain the Point2dC vector for it
		first_eyes = img_eyes.Copy();
		first_img = img.Copy();
		prev_snaxels = ArrayToSArray_B(pf_pts);
		SArray1dC<Point2dC> respts = AffineTransformPoints(prev_snaxels,aff_tr.Data2());
		DrawCrosses(img,respts);
	}
	

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
