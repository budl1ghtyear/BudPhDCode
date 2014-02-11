//File: HSPropJMeas

//Process: Uses particle filtering (HS Optical flow propagation, J based measurement)
//Author: Bud Goswami
//Date: 23.03.09
#include "Ravl/Array1d.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Image/DeinterlaceStream.hh"
#include "Ravl/Ellipse2d.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/DP/SequenceIO.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/DrawFrame.hh"
#include "Ravl/Image/OpenCVConvert.hh"
#include "Ravl/Image/ImageRectangle.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/IO.hh"
#include "Ravl/MeanCovariance.hh"
#include "Ravl/Option.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/String.hh"
#include "Ravl/SumsNd2.hh"
#include "opencv/highgui.h"
#include "opencv/cv.h"
//BUD SW
#include "ColourConvert.hh"//Contains colour conversion functions
#include "MouthRegion.hh"	//Contains mouth-region identification Functions 
#include "CostFunctions.hh"//Contains functions related to cluster grouping
#include "CascadedMCD.hh"
#include "EigenStateVector.hh"
#include "Particle.hh"
#include "SystematicResampling.hh"
#include "SIR.hh"
#include "TrackingFunctions.hh"
#include "SOARPropagationModel.hh"
#include "LLMeasurementModel.hh"
#include "EigenLoadParams.hh"
#include "TrackingFunctions.hh"
using namespace RavlImageN;
using namespace OmniN;
using namespace RavlConstN;
using namespace RavlN;
RandomGaussC BasePropagationModelC::rnd;

int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	DirectoryC qimg_dir = opt.String("i","TestSet/","Input Query Image Directory");
	FilenameC cascade_name = opt.String("c","/usr/share/opencv/haarcascades/haarcascade_frontalface_alt.xml","Path for the Haar Classifier");
	IntT segtype = opt.Int("t",1,"Type of Cluster Identification 0 - DensityCG, 1 - JCC");
	UIntT search = opt.Int("s",20,"Length of search window in pixels");
	UIntT noise = opt.Int("n",5,"Length of noise perturbtaion");
	RealT h = opt.Real("h",0.75,"h-value for MCD algorithm"); 
	IntT numpts = opt.Int("p",200,"Number of Particles");
	FilenameC mean_f = opt.String("mean","/vol/vssp/lip-tracking/LipShapeModels/TrainingData/OmniMean.txt","Path for the Model Mean");
	FilenameC evals_f = opt.String("evals","/vol/vssp/lip-tracking/LipShapeModels/TrainingData/OmniEVals.txt","Path for the Model EigenValues");
	FilenameC evect_f = opt.String("evect","/vol/vssp/lip-tracking/LipShapeModels/TrainingData/OmniEVect.txt","Path for the Model EigenVectors");
	FilenameC a1_f = opt.String("a1","/home/cvsspraid/cvsspraid2/cvssppg/eep1bg/workspace/PFTrackers/013A1Eigen.txt","Path for the A1 AR Parameters");
	FilenameC a2_f = opt.String("a2","/home/cvsspraid/cvsspraid2/cvssppg/eep1bg/workspace/PFTrackers/013A2Eigen.txt","Path for the A2 AR Parameters");
	FilenameC d_f = opt.String("d","/home/cvsspraid/cvsspraid2/cvssppg/eep1bg/workspace/PFTrackers/013DEigen.txt","Path for the D AR Parameters");
	opt.Compulsory("i");
	opt.Check();
	
	DListC<RealT> scorelist;	
	DListC<StringC> qfile = qimg_dir.FiltList("*.ppm");
	/////////////////////////////////////////////////
	//cout<<"LOAD THE MODEL PARAMETERS AND AR PARAMETERS"<<endl;
	MatrixC mean = GetMatrixNoSize(mean_f, 22,1);
	MatrixC evals = GetMatrixNoSize(evals_f, 22,1);
	MatrixC evect = GetMatrixNoSize(evect_f, 22,22);
	//cout<<"AR Parameters"<<endl;
	MatrixC a1 = GetMatrixWithSize(a1_f);
	MatrixC a2 = GetMatrixWithSize(a2_f);
	MatrixC d = GetMatrixWithSize(d_f);

	/////////////////////////////////////////////////
	cout<<"PERFORM BOOTSTRAP PROCESS ON FIRST IMAGE"<<endl;
	FilenameC first_img = qimg_dir + qfile.First();
	ImageC<RealRGBValueC> first_src;
	if(!Load(first_img, first_src)) cerr<<"Could not load the file "<<first_img<<endl;
	//Obtain face rec
	IndexC window = 50;
	ImageRectangleC botfacehalf = GetMouthRegion(first_img, cascade_name, window);
	Tuple2C<Point2dC,Point2dC> eyes_first = GetEyeLoc(first_img);
	DListC<Index2dC> lip_pix;	
	//For tracking, we are just going to use Normalised RG Image with CMCD
	ImageC<VectorC> m_img(GetNormRGVals(first_img),botfacehalf);
  	CascadedMCD mcd(m_img,h);
  	DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > finalpop = mcd.Apply();
  	lip_pix = ClusterIdent(finalpop, m_img, segtype);
  	ImageC<UIntT> lip_bin = GetBinaryLipImage(lip_pix,botfacehalf);
  	if(!Save("@X:LB Image",lip_bin)) cerr<<"Could not save Lip Binary Image"<<endl;
  	RealRGBValueC color(10,200,10); RealRGBValueC color2(10,10,200);
  	MeanCovarianceC lp_mc = GetMC(lip_pix,m_img);
  	MeanCovarianceC sk_mc = finalpop.First().Data1();
  	ImageC<RealT> j_img = ComputeJ(m_img, sk_mc, lp_mc);
  	Array1dC<Point2dC> lip_pt = GetLipContourPoints(lip_pix, 14, lip_bin, search, j_img);
  	cout<<"In Main Loop Now"<<endl;
  	for(Array1dIterC<Point2dC> it(lip_pt); it; it++)
  	{
  		if(first_src.Frame().Contains((*it)))
  			DrawCross(first_src,color,(*it),4);
  	} 

  	//Now we have the lip points, generate a B-Spline representation
  	//Wrap lip points
  	Array1dC<Point2dC> f_lip(1); f_lip[0] = lip_pt[0].Copy(); lip_pt.Append(f_lip);
  	UIntT order = 3; UIntT ncp = 11;
 	BSplineC bspl(order,ncp, BSplineC::UOPEN);
	Array1dC<Point2dC> cpts = bspl.CalculateControlPoints(lip_pt, BSplineC::CHORDLENGTH, BSplineC::UOPEN,order, ncp);	 	
	Array1dC<Point2dC> spl_pts = bspl.RenderCurve(30);
	for(Array1dIterC<Point2dC> it(spl_pts); it; it++)
	{
		if(first_src.Frame().Contains((*it)))
			DrawCross(first_src,color2,(*it),3);
	}
  	if(!Save("@X:Ellipse Image",first_src)) cerr<<"Could not save Lip Binary Image"<<endl;	

 	/////////////////////////////////////////////////
	cout<<"LOAD STREAMS"<<endl;
	DPIPortC<ImageC<RealRGBValueC> > inputstream;
	if(!OpenISequence(inputstream, (StringC)first_img))
	{
		cerr<<"Failed to open input sequence."<< qimg_dir<<" \n " ;
		return 1;
	}
	DeinterlaceStreamC<RealRGBValueC> dinputstream(inputstream);
	ImageC<RealRGBValueC> prev_img = dinputstream.Get();
	UIntT nth = 0;

  	/////////////////////////////////////////////////
  	cout<<"PARTICLE FILTER INSTANTIATION METHOD"<<endl;
  	//Get Head Affine Transform
  	ImageC<RealT> real_prev = ConvertRGBToRealImage(prev_img);
  	Tuple2C<ImageC<RealT>,Affine2dC> h_aff = ConvertImages(real_prev,eyes_first.Data1(),eyes_first.Data2());
	EigenStateVectorC initpts(cpts, h_aff.Data2(),mean, evect,noise);
	BaseStateVectorC *bsv = &initpts;
	EigenStateVectorC thetanmintwo = initpts;
	//cout<<"Test 1 = "<<cpts<<endl;
	//cout<<"Test 2 = "<<thetanmintwo.GetObservationVector()<<endl;
	//cout<<"Test 3 = "<<initpts.GetObservationVector()<<endl;
	//Instantiate Particle Filter Variables
	RealT initwt = 1.0 / (RealT) numpts;
	ParticleC initparticle(initpts,initwt);
	//ParticleC initparticle(,initwt);
	Array1dC<ParticleC> ptcol(numpts);
	ptcol.Fill(initparticle);
	//TRACKING PROCESS
	while(!dinputstream.IsGetEOS())//is end of stream?
	{
		ImageC<RealRGBValueC> img = dinputstream.Get();
		//Obtain affine transform
		FilenameC cur_file =  qimg_dir + qfile.Nth(nth++);
		Tuple2C<Point2dC,Point2dC> eyes;
		try
		{
			eyes = GetEyeLoc(cur_file);
		}
		catch(OmniN::ColossusExceptionNoFaceFoundC er)
  		{
    		eyes = eyes_first;		
  		}
  		cout<<"Got Eyes"<<endl;
  		//Get Affine Transform for head
		real_prev = ConvertRGBToRealImage(img);
  		h_aff = ConvertImages(real_prev,eyes.Data1(),eyes.Data2()); 		
		//INSERT TRACKING CODE IN HERE
		//Instantiate Propagation Model
		cout<<"Propagation"<<endl;
		SOARPropagationModelC ofprop(a1,a2,d,thetanmintwo.GetStateVector(),10);
		cout<<"Measurement"<<endl;
		ImageC<VectorC> curm_img(GetNormRGVals(img),botfacehalf.Expand(20));
		ImageC<RealT> j_img = ComputeJ(curm_img,sk_mc, lp_mc );
		LLMeasurementModelC ssdmeas(j_img,lp_mc,sk_mc,search);
		//Instantiate Resampling Model
		SystematicResamplingC sysres;
		//Check if particle collection is actually changing
		cout<<"PARTICLE COLLECTION "<<ptcol.Size()<<endl;
		for(UIntT i = 0; i < 5; i++)
		{
			cout<<ptcol[i].GetStateVector();
		}
		SIR pfilter(numpts, ofprop, ssdmeas, sysres,ptcol);
		cout<<"Result"<<endl;
		//BaseStateVectorC result = pfilter.Apply().Copy();
		//BaseStateVectorC *pf_res = pfilter.Apply();
		//Array1dC<Point2dC> pts = pf_res->GetObservationVector();
		Array1dC<Point2dC> pts = (pfilter.Apply()).GetObservationVector();
		cout<<"Vec - "<<pts<<endl;
		//Array1dC<Point2dC> pts = ConvertVecToPointArray(result.GetX());
		//manual affine transform
		cout<<"B-Spline Control Points Found were - "<<pts<<endl;
		//We now have B-Spline Control Points Returned To Us
		bspl.SetControlPoints(pts);
		Array1dC<Point2dC> spl_pts2 = bspl.RenderCurve(30);
		ImageC<RealRGBValueC> to_draw(img.Copy());
		for(Array1dIterC<Point2dC> it(spl_pts2); it; it++)
		{
			Point2dC pt((*it).Row(), (*it).Col());
	  		if(to_draw.Frame().Contains(pt))
  				DrawCross(to_draw,color,pt,4);		
		}
		//Now we have a result for the Particle Filter Estimate, just obtain the Point2dC vector for it
		if(!Save("@X: Str Image",to_draw)) cerr<<"Could not save stream image"<<endl;
		EigenStateVectorC dummy(pts, h_aff.Data2(),mean, evect, 5);
		thetanmintwo = dummy;
		//IntT s;
		//cin>>s;
	}	
	return 0;
}

