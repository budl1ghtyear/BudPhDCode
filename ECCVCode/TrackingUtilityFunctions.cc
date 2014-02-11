#include "TrackingUtilityFunctions.hh"
/*
Tuple2C<BasisSplineC,Affine2dC> GetLipPoints(FilenameC const &lipfile)
{
	//This method serves to return the co-ordinates of the lip given a filename
	ImageC<RealRGBValueC> img;
	if(!Load(lipfile, img)) cerr<<"Could not load input image"<<endl;	
	return GetLipPoints(img);
}

Tuple2C<BasisSplineC,Affine2dC> GetLipPoints(ImageC<RealRGBValueC> const &lipfile)
{
	ImageC<RealRGBValueC> img(lipfile.Copy());
	if(!Save("@X:First Image",img)) cerr<<"Could not display the results of the bootstrap process"<<endl;
	LipSegmentationC lipseg(lipfile);
	//"ClusteringType - values allowed - CMCD=0,KSMCD=1,KM=2,FCM=3, default = CMCD"
	IntT cluster_type = 0;
	//"Pixel Type - RGB=0, NormRG=1, HSV=2, YUV=3, CIE=4, PHue=5 , default = RGB=0"
	IntT pixel_type = 4;
	//"Pixel labelling type - JNormal=0, JSkin=1, default =  JNormal=0 "
	IntT label_type = 0;
	//"Region Identification type - CC=0,CCentral=1, default =CC=0 "
	IntT region_type = 1;
	RealT hvalue = 0.8;
	UIntT numclust = 3;
	//At this stage we have an identified region image containing the lip region
	Tuple3C<ImageC<UIntT>, MeanCovarianceC, MeanCovarianceC> result = lipseg.Apply((PixelType)pixel_type,(ClusteringType)cluster_type,(LabellingType) label_type,(RegionIdentificationType) region_type,hvalue, numclust);
	cout<<"BootStrap Im Frame = "<<result.Data1().Frame()<<endl;
	if(!Save("@X:Identified Lip Region",result.Data1())) cerr<<"Could not display the results of the bootstrap process"<<endl;
	sleep(5);
	ImageC<UIntT> lip_img = MakeFullSizeImage(result.Data1(), img.Frame());
	LipExtractorC lext(lip_img, result.Data1().Frame(), 15);
	SArray1dC<Point2dC> lippts = lext.GetLipBoundaryPoints();
	cout<<"Obtained Lip Boundary Points "<<lippts<<endl;
	BasisSplineC bspl(11,2,30);
	bspl.FitCurve(lippts);
	//Draw the result on the provided image:
	SArray1dC<Point2dC> curvepts = bspl.GetCurvePoints();
	RealRGBValueC cross(250,10,10);
	for(SArray1dIterC<Point2dC> it(curvepts); it; it++)
	{
		DrawCross(img,cross,Index2dC((*it).Row(),(*it).Col()));
	}
	if(!Save("@X:Generated Spline Points For First Image",img)) cerr<<"Could not display the results of the bootstrap process"<<endl;
	Tuple2C<BasisSplineC,Affine2dC> res(bspl,lipseg.GetAffineTransform().Inverse());
	return res;	
}

SArray1dC<Point2dC> AffineTransformPoints(const SArray1dC<Point2dC> &pts, const Affine2dC &aff)
{
	SArray1dC<Point2dC> res(pts.Size()); res.Fill(Point2dC(0.0,0.0));
	for(SArray1dIter2C<Point2dC, Point2dC> it(pts,res); it; it++)
	{
		Vector2dC vec(it.Data1().Row(),it.Data1().Col());
		vec = aff * vec;
		Point2dC pt(vec.Row(),vec.Col());
		it.Data2() = pt.Copy();
	}
	return res;
}
* */

ImageC<RealT> ConvertRGBToRealImage(const ImageC<RealRGBValueC> &im)
{
	//Conversion uses:
	//0.2989 * R + 0.5870 * G + 0.1140 * B 
	ImageC<RealT> res(im.Frame(),0.0);
	for(Array2dIter2C<RealRGBValueC,RealT> it(im, res); it; it++)
	{
		it.Data2() = (0.2989 * it.Data1().Red()) + (0.5870 * it.Data1().Green()) + (0.1140 * it.Data1().Blue());
	}
	return res;
}


ImageC<RealT> GetSusanImage(ImageC<RealT> const &img) //PLEASE BEAR IN MIND THAT RUNNING THIS METHOD EXPRESSLY REQUIRES A COMPILED INSTANCE OF THE SUSAN PROGRAM TO BE AVAILABLE
{
	FilenameC infile = "/dev/shm/inimg.pgm";
	if(!Save(infile,img)) cerr<<"Could not save the input image to susan"<<endl;
	FilenameC outfile = "/dev/shm/outimg.pgm";
	ImageC<RealT> result;
	//Now run the system command to perform the susan operation
	StringC cmd = "~/Code/Susan/susan "+infile+" "+outfile+" -se";
	if(system(cmd))
	{
		cerr<<"Could not execute the system command : "<<cmd<<endl;
		exit(1);
	}
	else
	{
		if(!Load(outfile,result))
			cerr<<"Could not load the results file"<<endl;
	}
	//Now remove the written images
	if(((system("rm "+infile))&&(system("rm "+outfile))))
	{
		cerr<<"Could not remove the result of the susan image manipulation"<<endl;
		exit(1);
	}
	//Normalise the resulting image
	return NormaliseImage(result);
}
/*
Tuple2C<ImageC<RealT>, ImageC<Vector2dC> > GetMeasurementImages(ImageC<RealRGBValueC> const &im1, Tuple2C<Point2dC,Point2dC> const &eyes1,ImageC<RealRGBValueC> const &im2, Tuple2C<Point2dC,Point2dC> const &eyes2)
{
	//Image One is the current image
	//Image Two is the previous image
	Tuple2C<ImageC<RealRGBValueC>,Affine2dC> aff_tr1 = AffNormWRTEyes(im1,eyes1.Data1(),eyes1.Data2()); //create the first affine transformed image
	Tuple2C<ImageC<RealRGBValueC>,Affine2dC> aff_tr2 = AffNormWRTEyes(im2,eyes2.Data1(),eyes2.Data2());
	ImageRectangleC botfacehalf(100,aff_tr1.Data1().BRow(),20,110); //the necessary image rectangle
	ImageC<RealRGBValueC> rim1(aff_tr1.Data1(),botfacehalf); ImageC<RealT> rim1real = ConvertRGBToRealImage(rim1);
	ImageC<RealRGBValueC> rim2(aff_tr2.Data1(),botfacehalf); ImageC<RealT> rim2real = ConvertRGBToRealImage(rim2);
	//Perform Susan Based Edge Detection
	ImageC<RealT> susanim = GetSusanImage(rim1real);
	//Now perform Optical Flow Measurement
	LMSOpticFlowC of; PairC<ImageC<RealT> > ofims(rim1real,rim2real);
	ImageC<Vector2dC> ofim = of.Estimate(ofims);
	return Tuple2C<ImageC<RealT>, ImageC<Vector2dC> >(susanim,ofim);	
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
*/
void DisplayImage(const SArray1dC<Point2dC> &pts, const UIntT &size)
{
	ImageC<RealRGBValueC> img(IndexRange2dC(0,700,0,700));
	img.Fill(RealRGBValueC(0.0,0.0,0.0));
	//Generate Random Colour
	RealRGBValueC rgb(255,10,10);
	for(SArray1dIterC<Point2dC> it(pts); it; it++)
	{
		DrawCross(img,rgb,Index2dC((*it).Row(),(*it).Col()),size);
	}
	if(!Save("@XA:Output Image",img)) cerr<<"Could not save the resulting output image"<<endl;
}

ImageRectangleC ObtainMouthLocalisationWindow(SArray1dC<Point2dC> const &pts, UIntT const &src)
{
	//Obtain the various parameters for trow, brow, rcol and lcol
	RealT trow(1000.),brow(0.),lcol(1000.),rcol(0.);
	for(SArray1dIterC<Point2dC> it(pts); it; it++)
	{
		if(it->Row() < trow)
			trow = it->Row();
		if(it->Row() > brow)
			brow = it->Row();
		if(it->Col() < lcol)
			lcol = it->Col();
		if(it->Col() > rcol)
			rcol = it->Col();		
	}
	ImageRectangleC imrec(trow,brow,lcol,rcol);
	return imrec.Expand(IndexC(src));
}

