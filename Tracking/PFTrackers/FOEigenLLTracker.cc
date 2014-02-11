//LipClusteringC test program

//Just to test the Lip Clustering Test Class
#include "Ravl/Affine2d.hh"
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/SArray1d.hh"
#include "Ravl/SArray1dIter3.hh"
#include "Ravl/Option.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/DP/SequenceIO.hh"
#include "Ravl/Image/DeinterlaceStream.hh"
#include "Ravl/String.hh"
#include "LipSegmentation.hh"
#include "ColourSpaceImage.hh"
#include "XMLToGTConv.hh"
#include "BANCAGT.hh"
#include "Ravl/Boundary.hh"
#include "Ravl/Image/MorphClose.hh"
#include "Ravl/Image/Erode.hh"
#include "Ravl/Ellipse2d.hh"
#include "Ravl/Image/DrawEllipse.hh"
#include "Ravl/Array2dIter3.hh"
#include "Ravl/OS/Date.hh"
#include "Ravl/Point2d.hh"
#include "LipExtractorC.hh"
//#include "BSplineSolver.hh"
#include "BSplineC.hh"
//#include "PointStateVector.hh"
#include "Particle.hh"
#include "SystematicResampling.hh"
#include "SIR.hh"
#include "TrackingFunctions.hh"
#include "FOPropagationModel.hh"
#include "LLMeasurementModel.hh"
#include "TrackingFunctions.hh"
#include "EAStateVector.hh"
//#define NUM_POINTS_EACH_SIDE 6
using namespace RavlN;
using namespace RavlConstN;
using namespace RavlImageN;
RandomGaussC BasePropagationModelC::rnd;
RealT ComputeSegmentationQuality(const ImageC<UIntT> &im, const ImageC<UIntT> &gt);
//Method to compute the image quality based on binary ground truth and segmented lip image
Array1dC<Point2dC> AffineTransformPoints(const Array1dC<Point2dC> &pts, const Affine2dC &aff);
int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	DirectoryC img_dir = opt.String("i","Input/","QUERY INPUT DIRECTORY");
	StringC file_ext = opt.String("e","ppm","QUERY INPUT FILE EXTENSION");
	UIntT search = opt.Int("s",20,"Length of search window in pixels in affine space");
	UIntT noise = opt.Int("n",2,"Length of noise perturbtaion in pixels");
	IntT numpts = opt.Int("p",200,"Number of Particles");
	RealT hvalue = opt.Real("h",0.75,"Clustering Confidence Value");
	opt.Compulsory("i");
	opt.Check();

	//SEGMENTATION DEFAULT VALUES:
	IntT cluster_type = 0; IntT pixel_type = 0; IntT label_type = 0; IntT region_type = 1; UIntT ncl = 3;
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
	
	IStreamC inMean(StringC("/vol/vssp/lip-tracking/LipShapeModels/TrainingData/OmniMean.txt"));
	RealT meanValue;
	DArray1dC<RealT> meanValues;
	while (!inMean.eof())
	{	inMean>>meanValue;
		meanValues.Append(meanValue);
	}
	Array1dC<Point2dC> means(meanValues.Size()/2);
	IndexC i=0;
	for (i=0;i<means.Size();i++)
	{
		means[i]= Point2dC(meanValues[i+means.Size()],meanValues[i]);
	}
	inMean.Close();
	cout << means << endl;
	
	IStreamC inEVect(StringC("/vol/vssp/lip-tracking/LipShapeModels/TrainingData/OmniEVect.txt"));
	RealT eVectValue;
	MatrixC eVects(2*means.Size(),2*means.Size());
	cout <<" Line 134 "<< eVects.Size() << endl;
	i=0;
	UIntT pager = 2*means.Size() ;
	while (!inEVect.eof())
	{	
		RealT temp;
		inEVect>>temp;
		if (i<(2*means.Size())*(2*means.Size()))
		{
			eVects[i/pager][i % pager]=temp;
			//if (i % pager == 0) cout << i << endl;
			i++;
		}
	}
	inEVect.Close();
	//cout << eVects << endl;

	
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
	LipSegmentationC lipseg(first_img);
	Tuple3C<ImageC<UIntT>, MeanCovarianceC, MeanCovarianceC> result = lipseg.Apply((PixelType)pixel_type,(ClusteringType)cluster_type,(LabellingType) label_type,(RegionIdentificationType) region_type,hvalue,ncl);
	cout<<"Mean Covariance Parameters "<<result.Data2()<<"\n"<<result.Data3()<<endl;
	//Make full size Binary Lip Image
	ImageC<UIntT> lip_img(first_img.Frame(),0);
	DListC<Index2dC> lip_ind;
	for(Array2dIterC<UIntT> it(result.Data1()); it; it++)
	{
		Index2dC ind = it.Index();
		if((*it) == 255)
		{
			lip_img[ind] = 255;
			lip_ind.Append(ind.Copy());
		}
	}
	LipExtractorC lext(lip_img,result.Data1().Frame(),20);
	Array1dC<Point2dC> lip_points(lext.GetLipBoundaryPointsArray());
	//cout<<"Lip Points "<<lip_points<<endl; 
	BSplineC bspl(4, 11, BSplineC::UOPEN, false); //REMEMBER - CUBIC B-SPLINE WITH 11 CONTROL POINTS
	Array1dC<Point2dC> spl_cp = bspl.CalculateControlPoints(lip_points,BSplineC::UOPEN,BSplineC::UNIFORM); //UNIFORM OPEN KV, WITH UNIFORM PARSEL
	// THIS PROCESS IS STANDARD FOR ALL PARTICLE FILTERS UP TILL LINE 62opt.String("i","Input/","QUERY INPUT DIRECTORY");
	// NOW WE HAVE CONTROL POINTS, INSTANTIATE THE PARTICLE FILTER 
	// REMEMBER TO DO THIS IN AFFINE SPACE, SO TRANSFORM TO AFFINE SPACE
	Array1dC<Point2dC> aff_spl_cp = AffineTransformPoints(spl_cp,lipseg.GetAffineTransform().Inverse());
	//PointStateVectorC initpts(aff_spl_cp);
	//BSplineStateVectorC initpts(aff_spl_cp,bspl.GetKnotVector(),20,4); // Instantiate new B-Spline State Vector Class
	EAStateVectorC initpts(means,eVects,lipseg.GetAffineTransform().Inverse(),aff_spl_cp,bspl.GetKnotVector());
	RealT initwt = 1.0 / (RealT) numpts;
	ParticleC initparticle(initpts,initwt);
	Array1dC<ParticleC> ptcol(numpts);
	ptcol.Fill(initparticle);	
	// NOW BASIC PARTICLE COLLECTION HAS BEEN INSTANTIATED
	UIntT nth = 1;
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
		//NO AFFINE TRANSFORM THE TRACKING IMAGE
		Tuple2C<ImageC<RealRGBValueC>,Affine2dC> aff_tr = AffNormWRTEyes(img,img_eyes.Data1(),img_eyes.Data2());
		FOPropagationModelC ofprop(noise); //PROPAGATION MODEL INSTANTIATED
		ImageRectangleC botfacehalf(100,aff_tr.Data1().BRow(),20,110);
		ImageC<RealRGBValueC> curm_img(aff_tr.Data1(),botfacehalf);
		cout<<"Compute J Matrix "<<endl;
		ImageC<RealT> j_img = ComputeJ(curm_img,result.Data2(),result.Data3());
		LLMeasurementModelC ssdmeas(j_img,result.Data3(),result.Data2(),search);
		SystematicResamplingC sysres;
		SIR pfilter(numpts, ofprop, ssdmeas, sysres,ptcol);
		cout<<"Starting Particle Filter here"<<endl;
		Array1dC<Point2dC> pf_pts = pfilter.Apply().GetObservationVector(); //GET THE APPLY METHOD 
		Array1dC<Point2dC> img_pts = AffineTransformPoints(pf_pts,aff_tr.Data2());
		BSplineC spl(bspl.GetOrder(),bspl.GetKnotVector(),img_pts); //DRAW THE SPLINE USING OUR INTIALISED VECTOR
		ImageC<RealRGBValueC> to_draw(img.Copy());
		for(Array1dIterC<Point2dC> it(img_pts); it; it++)
		{
			Point2dC pt((*it).Row(), (*it).Col());
	  		if(to_draw.Frame().Contains(pt))
  				DrawCross(to_draw,color,pt,4);		
		}
		//Now we have a result for the Particle Filter Estimate, just obtain the Point2dC vector for it
		if(!Save("@X: Str Image",to_draw)) cerr<<"Could not save stream image"<<endl;
		first_eyes = img_eyes.Copy();
	}
	return 0;
}

Array1dC<Point2dC> AffineTransformPoints(const Array1dC<Point2dC> &pts, const Affine2dC &aff)
{
	Array1dC<Point2dC> res(pts.Size()); res.Fill(Point2dC(0.0,0.0));
	for(Array1dIter2C<Point2dC, Point2dC> it(pts,res); it; it++)
	{
		Vector2dC vec(it.Data1().Row(),it.Data1().Col());
		vec = aff * vec;
		Point2dC pt(vec.Row(),vec.Col());
		it.Data2() = pt.Copy();
	}
	return res;
}

