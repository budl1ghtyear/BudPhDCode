//LipClusteringC test program
//LipClusteringC test program

//Just to test the Lip Clustering Test Class
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/Option.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/String.hh"
#include "Ravl/Tuple2.hh"
#include "LipSegmentation.hh"
#include "ColourSpaceImage.hh"
#include "XMLToGTConv.hh"
#include "BANCAGT.hh"
#include "BSplineC.hh"
#include "LipExtractorC.hh"
#include "Ravl/Image/MorphClose.hh"
#include "Ravl/Image/Erode.hh"
#include "Ravl/Image/MorphOpen.hh"
#include "CostFunctions.hh"
using namespace RavlN;
using namespace RavlImageN;
using namespace OmniN;
//RealT ComputeSegmentationQuality(const ImageC<UIntT> &im, const ImageC<UIntT> &gt);
//Method to compute the image quality based on binary ground truth and segmented lip image
ImageC<UIntT> MakeFullSizeImage(const ImageC<UIntT> &r_img, const ImageRectangleC &frame);
ImageC<UIntT> RegionIdentificationEllipseMassDensity(const ImageC<UIntT> &m_img);
//bool IsInsideEllipse(const Ellipse2dC &ell, const Index2dC &ind);
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
	if(!Save("@X:Original Image",img)) cerr<<"Could not show input image"<<endl;
	LipSegmentationC lipseg(img_name);
	//Colour Space Conversion
	ColourSpaceImage col_img = lipseg.ColourSpaceConversion((PixelType)pixel_type);
	//Affine Mouth Projection
	cout<<"Colour space conversion performed"<<endl;
	Tuple4C<Affine2dC, ImageRectangleC, ImageRectangleC,ImageC<VectorC> > amri = lipseg.AffineMouthRegionDetection(col_img);
	cout<<"Affine Mouth Region Detection Performed"<<endl;
	//Clustering
	Tuple2C<MeanCovarianceC,MeanCovarianceC> col_trends = lipseg.LipClusteringStep((ClusteringType)cluster_type,ImageC<VectorC>(amri.Data4(),amri.Data2()),hvalue,numclust);
	cout<<"Lip Clustering Step Performed"<<endl;
	/*
	 * //Draw Clustered Image
	ImageC<UIntT> clust_img(img.Frame(),0); UIntT label = 0;
	for(DLIterC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > it(col_trends); it; it++)
	{
		for(DLIterC<Tuple2C<VectorC,Index2dC> > it2((*it).Data2()); it2; it2++)
		{
			clust_img[(*it2).Data2()] = label;
		}
		label++;
	}
	SegmentationC seg(clust_img,label); 
	if(!Save("@X:Clustered Image",seg.RandomImage())) cerr<<"Could not save the clustered image"<<endl;
	*/
	//Now we should have the cluster data available to us, we can change the pixel labelling types and region identification types as we see fit
	//Perform scoring with l0 and r0
	ImageC<UIntT> ppl_img = lipseg.ProbabilisticPixelLabellingStep(ImageC<VectorC>(col_img.GetImage(),amri.Data3()) , col_trends.Data1(), col_trends.Data2(), (LabellingType)0);
	ImageC<UIntT> rgn_img = lipseg.RegionIdentificationStep(ppl_img,(RegionIdentificationType)0);
	ImageC<UIntT> lip_img = MakeFullSizeImage(rgn_img, img.Frame());
	if(!Save("@X:PPL Image",ppl_img)) cerr<<"Could not save the clustered image"<<endl;
	if(!Save("@X:RGI Image",rgn_img)) cerr<<"Could not save the clustered image"<<endl;
	if(!Save("@X:BINARY Image",lip_img)) cerr<<"Could not save the clustered image"<<endl;
	Tuple2C<MeanCovarianceC,MeanCovarianceC> col = lipseg.GetColorTrends();
	ImageC<UIntT> test_img(ppl_img.Frame()); test_img.Fill(255);
	ImageC<VectorC> test_col_img(col_img.GetImage(), ppl_img.Frame());
	for (Array2dIter2C<VectorC, UIntT> it(test_col_img, test_img); it; it++)
	{
		if(col.Data2().MahalanobisDistance(it.Data1()) < 6.0)
			it.Data2() = 0;
	}
	if(!Save("@X:NON-SKIN Image",test_img)) cerr<<"Could not save the clustered image"<<endl;
	//Maybe perform morphology on PPL image to remove salt and pepper noise
	ImageC<UIntT> kernel(IndexRangeC(-1,1),IndexRangeC(-1,1));
	kernel.Fill(1);
	ImageC<UIntT> opened_img,salt_img;
	//Remove salt and pepper noise - opening followed by closing
	MorphBinaryOpen2d(ppl_img,kernel,opened_img,(UIntT)1);
	MorphBinaryClose2d(opened_img,kernel,salt_img,(UIntT)1);
	if(!Save("@X:opened Image",opened_img)) cerr<<"Could not save the clustered image"<<endl;
	if(!Save("@X:opened,closed Image",salt_img)) cerr<<"Could not save the clustered image"<<endl;
	//Perform connected components analysis on the resulting image
	/*ConnectedComponentsC<UIntT> connected ; 
	Tuple2C<ImageC<UIntT>, UIntT> result = connected.Apply(salt_img); 
	SegmentationC seg_map(result);
	if(!Save("@X:SegMap Image",seg_map.RandomTaintImage())) cerr<<"Could not save the clustered image"<<endl;
	//get the HSV equivalent image
	ImageC<RealHSVValueC> hsv_img = RealRGBImageCT2RealHSVImageCT(ImageC<RealRGBValueC>(img,ppl_img.Frame()));
	if(!Save("@X:HSV Image",hsv_img)) cerr<<"Could not save the clustered image"<<endl;*/
	ImageC<UIntT> res_ell_img  = RegionIdentificationEllipseMassDensity(salt_img);
	if(!Save("@X:ELLMI Image",res_ell_img)) cerr<<"Could not save the clustered image"<<endl;

	/*
	//B-Spline Stuff
	LipExtractorC lext(lip_img,rgn_img.Frame(),20);
	Array1dC<Point2dC> lip_points(lext.GetLipBoundaryPointsArray());
	cout<<"Lip Points "<<lip_points<<endl; 
	//BSplineSolver bspl(lip_points, 4, 11, BSplineSolver::AVERAGING, BSplineSolver::CHORDLENGTH);
	BSplineC bspl(4, 11, BSplineC::UOPEN, false);
	Array1dC<Point2dC> spl_cp = bspl.CalculateControlPoints(lip_points,BSplineC::UOPEN,BSplineC::UNIFORM);
	//BSplineC bspl(4, 11, BSplineC::UPCLOSED, false);
	//Array1dC<Point2dC> spl_cp = bspl.CalculateControlPoints(lip_points,BSplineC::UPCLOSED,BSplineC::UNIFORM);
	cout<< "B-Spline Control Points Are - \n"<< spl_cp << endl;
	//Array1dC<Point2dC> spl_pts = bspl.RenderCurve(500);
	Array1dC<Point2dC> spl_pts = bspl.RenderCurve(100);
	ImageC<RealRGBValueC> spl_img(img.Copy());
	RealRGBValueC cp_col(20,255,20);RealRGBValueC pt_col(255,20,20); RealRGBValueC lip_col(20,20,20);
	
	for(Array1dIterC<Point2dC> it(spl_pts); it; it++)
	{
		DrawCross(spl_img,pt_col,Index2dC((*it).Row(),(*it).Col()),1);
	}
	
	for(Array1dIterC<Point2dC> it(spl_cp); it; it++)
	{
		DrawCross(spl_img,cp_col,Index2dC((*it).Row(),(*it).Col()),3);
	}
	for(Array1dIterC<Point2dC> it(lip_points); it; it++)
	{
		DrawCross(spl_img,lip_col,Index2dC((*it).Row(),(*it).Col()),4);
	}
	if(!Save("@X:Spline Lip Image",spl_img)) cerr<<"Could not save binary lip image"<<endl;	
	//Lip Extractor Label Images Stuff
	ImageC<UIntT> lext_img = lext.PerformMorphology().Copy();
	ImageC<RealRGBValueC> lext_img_draw = img.Copy(); RealRGBValueC contour_col(10,10,10);
 	for(Array2dIter2C<UIntT, RealRGBValueC> it(lext_img, lext_img_draw); it; it++)
	{
		if(it.Data1() == 1.0)
		{
			DrawCross(lext_img_draw,contour_col,it.Index(),1);
		}
	}
	if(!Save("@X:Drawn Contour Lip Image",lext_img_draw)) cerr<<"Could not save binary lip image"<<endl;	
	ImageC<RealRGBValueC> lip_reg_img(img, rgn_img.Frame());
	if(!Save("@X: Lip Region Image",lip_reg_img)) cerr<<"Could not show the lip region image"<<endl;
	//Get GT region
	/*
	FilenameC gt_img_name = gt_dir + img_name.BaseNameComponent() + ".pbm";
	ImageC<UIntT> gt_img;
	if(!Load(gt_img_name,gt_img)) cerr<<"Could not load Ground Truth Image"<<endl;
	ImageC<UIntT> draw_gt_img(gt_img,rgn_img.Frame());
	if(!Save("@X:GT Region",draw_gt_img)) cerr<<"Could not show GT image"<<endl;
	ImageC<RealRGBValueC> draw_cnt_img(lext_img_draw,rgn_img.Frame());
	if(!Save("@X:Contour Spec Region",draw_cnt_img)) cerr<<"Could not show GT image"<<endl;
	*/
	return 0;
}
/*
RealT ComputeSegmentationQuality(const ImageC<UIntT> &im, const ImageC<UIntT> &gt)
{
	RealT num = 0.0,den=0.0;
	for(Array2dIter2C<UIntT, UIntT> it(im,gt); it; it++)
	{
		if(it.Data1() == 255)
		{
			den++;
			if(it.Data2() == 255)
			{
				num++;
			}
		}
		else
		{
			 if(it.Data2() == 255)
			{
				den++;
			}
		}
	}
	return ((RealT)(num/den));
}
*/
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
ImageC<UIntT> RegionIdentificationEllipseMassDensity(const ImageC<UIntT> &m_img)
{
	//Perform Connected Components Analysis on this
	ConnectedComponentsC<UIntT> connected;
	Tuple2C<ImageC<UIntT>, UIntT> result = connected.Apply(m_img); 	
	//Tuple2C<ImageC<UIntT>, UIntT> result(m_img,2); 	cd Mode
	//Now we have a segmentation image with all zeros except for the lip and noise areas
	SegmentationC seg(result.Data1(),result.Data2());
	seg.RemoveSmallComponents(100);
	if(!Save("@X: Ellipse Segmentation Image",seg.RandomImage()))
		cerr<<"Failed to show segmentation image"<<endl;
	//Now we have a situation where we can compute the eccentricity of non-zero clusters
	Array1dC<DListC<Index2dC> > clusters(seg.Labels());
	for(Array2dIterC<UIntT> it(seg.SegMap()); it; it++)
	{
		clusters[*it].Append(it.Index().Copy());
	}
	RealT cost = 0.0;
	Index2dC centre = m_img.Frame().Center();
	DListC<Index2dC> lipind;ImageC<UIntT> res_img(m_img.Frame());res_img.Fill(0);
	for(Array1dIterC<DListC<Index2dC> > it(clusters); it; it++)
	{
		SArray1dC<Point2dC> elldat((*it).Size());
		IndexC i = 0;
		for(DLIterC<Index2dC> dt((*it)); dt; dt++)
		{
			elldat[i] = (Point2dC)(*dt);
			i++;
		}
		Ellipse2dC ell;
		FitEllipse(elldat,ell);
		//Weight = mass (points inside the ellipse)/ area of the cluster
		UIntT mass = 0;
		for(SArray1dIterC<Point2dC> s(elldat); s; s++)
		{
			Index2dC ellind((*s).Row(), (*s).Col());
			if(IsInsideEllipse(ell, ellind))
				mass++;
		}
		RealT f = 0.0;
		if((mass > 0.0)&&(it.Index() != 0))
			f = (RealT)mass / (RealT)(elldat.Size()*centre.SqrEuclidDistance((Index2dC)ell.Centre())) ; 
		if(f>cost)
		{
			cost = f;
			lipind = (*it).Copy(); 
		}			
	}		
	for(DLIterC<Index2dC> it(lipind); it; it++)
	{
		res_img[*it] = 255;
	}
	return res_img;		
}
/*
bool IsInsideEllipse(const Ellipse2dC &ell, const Index2dC &ind)
{
//~ I assume you also know the location of the ellipse's center. Call that (x0,y0).
//~ Let t be the counterclockwise angle the major axis makes with respect to the
//~ x-axis. Let a and b be the semi-major and semi-minor axes, respectively. If
//~ P = (x,y) is an arbitrary point then do this:
//~ 
 //~ X = (x-x0)*cos(t)+(y-y0)*sin(t); % Translate and rotate coords.
 //~ Y = -(x-x0)*sin(t)+(y-y0)*cos(t); % to align with ellipse
//~ 
//~ If
//~ 
 //~ X^2/a^2+Y^2/b^2
//~ 
//~ is less than 1, the point P lies inside the ellipse. If it equals 1, it is right on
//~ the ellipse. If it is greater than 1, P is outside.

	RealT ymyzer = (RealT)(ind.Row() - ell.Centre().Row());
	RealT xmxzer = (RealT)(ind.Col() - ell.Centre().Col());
	RealT maj,min,angle;
	Point2dC centre;
	if( ell.EllipseParameters(centre,maj,min,angle))
	{
		RealT X = (xmxzer*Cos(angle)) + (ymyzer*Sin(angle));
		RealT Y = (xmxzer*Sin(angle)) + (ymyzer*Cos(angle));
		maj /= 2.0;//Compute semi major and minor axes
		min /= 2.0;
		RealT value = (Sqr(X)/Sqr(maj))+(Sqr(Y)/Sqr(min));
		if(value <= 1.0)
			return true;
		else
			return false;
	}
	else
		return false;
}
*/

