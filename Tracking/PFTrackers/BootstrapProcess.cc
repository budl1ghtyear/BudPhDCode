//File: CascadedMCDScorer.cc
//Input: Query Image Directory, Output Image Directory, Type of Colour Space, Type of Cost Function
//Process: Uses Cascaded MCD to perform automatic lip segmentation
//Author: Bud Goswami
//Date: 23.03.09
#include "Ravl/Array1d.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/ConnectedComponents.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/DrawFrame.hh"
#include "Ravl/Image/OpenCVConvert.hh"
#include "Ravl/Image/ImageRectangle.hh"
#include "Ravl/Image/Segmentation.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/IO.hh"
#include "Ravl/PatternRec/DesignGaussianMixture.hh"
#include "Ravl/PatternRec/Sample.hh"
#include "Ravl/PatternRec/DistanceMahalanobis.hh"
#include "Ravl/PatternRec/Classifier.hh"
#include "Ravl/PatternRec/DesignKMeans.hh"
#include "Ravl/PatternRec/DesignFuzzyCMeansCluster.hh"
#include "Ravl/SumsNd2.hh"
#include "Ravl/Option.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/String.hh"
#include "opencv/highgui.h"
#include "opencv/cv.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <ctype.h>
#include "Omni/ColossusException.hh"
#include "Omni/FaceScan.hh"
#include "Omni/DetectFace.hh"
#include "Omni/FaceDetectionModel.hh"
#include "Omni/DetectedFace.hh"
#include "Omni/RawFace.hh"
#include "Omni/DisplayFace.hh"
//BUD SW
//#include "LoadingFunctions.hh"
#include "ColourConvert.hh"//Contains colour conversion functions
#include "MouthRegion.hh"	//Contains mouth-region identification Functions 
#include "CostFunctions.hh"//Contains functions related to cluster grouping
#include "CascadedMCD.hh"
#include "KSmirnovMCD.hh"
#include "BSplineC.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/DrawFrame.hh"
#include "Ravl/Image/DrawLine.hh"
//#include "TrackingFunctions.hh"
using namespace RavlImageN;
using namespace OmniN;
RandomGaussC BasePropagationModelC::rnd;
template <typename CL>
DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > MCDClustering(const ImageC<CL> &im, const UIntT &c_type, const RealT &hval);
template <typename T>
ImageC<RealT> ProbabilisticPixelLabellingStep(const ImageC<T> &im, const MeanCovarianceC &sk, const MeanCovarianceC &lp, const UIntT &type);
Array1dC<Point2dC> GetLipContourPoints(const DListC<Index2dC> &lp, const UIntT &res, const ImageC<UIntT> &img, const RealT &search, const ImageC<RealT> &jim);
SegmentationC CentralColumnRegions(const SegmentationC &seg);
ImageC<UIntT> PerformSpatialRegionIdentification(const ImageC<UIntT> &im, const ImageRectangleC &imrec, const UIntT &type);
ImageC<UIntT> GetLargestConnectedComponent(const ImageC<UIntT> &im, const ImageRectangleC &imrec);
int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	DirectoryC qimg_dir = opt.String("i","TestSet/","Input Query Image Directory with ppm images");
	DirectoryC gtqimg_dir = opt.String("g","GT/","Input Ground Truth Image Directory");
	DirectoryC outimg_dir = opt.String("o","Output/","Output Image Directory");
	FilenameC cascade_name = opt.String("c","/usr/share/opencv/haarcascades/haarcascade_frontalface_alt.xml","Path for the Haar Classifier");
	StringC extension = opt.String("e",".ppm","Input Image Extension");
	UIntT pix_type = opt.Int("p",0,"Pixel Output Type - 0=RealRGBValueC, 1=NormalisedRG, 2=RealHSVValueC, 3=RealYUVValueC, 4=CIELab, 5 = PseudoHue");
	UIntT clust_type = opt.Int("l",0,"Robust Clustering Type - 0=CascadedMCD, 1=MCD with KS, 2=K Means, 3= Fuzzy C Means");
	UIntT lab_type = opt.Int("pl",0,"Type of probabilistic pixel labelling- 0 = MskinSskinMlipSlip, 1 = MskinSskinMlipSskin");
	UIntT reg_type = opt.Int("r",0,"Type of spatial region identification- 0 = Largest Connected Component");
	RealT h = opt.Real("h",0.75,"Confidence Value for Methods");
	opt.Check();

	//RUN LOOP FOR DIRECTORY IMAGES
	DListC<RealT> scorelist;	
	DListC<StringC> qfile = qimg_dir.FiltList("*"+extension);
	for(DLIterC<StringC> dt(qfile); dt; dt++)
	{
		FilenameC qimg = qimg_dir+(*dt);		
		FilenameC oimg_seg = outimg_dir+"SegmentationResults/"+qimg.BaseNameComponent()+"SEG.ppm";
		FilenameC oimg_spl = outimg_dir+"SplineResults/"+qimg.BaseNameComponent()+"SPL.ppm";
		FilenameC gtimg = gtqimg_dir + qimg.BaseNameComponent() + ".pbm";
		//Load an example in RGB
		ImageC<RealRGBValueC> src;   
		if(!Load(qimg, src)) cerr<<"Loading RAVL Image Failed"<<endl;
		if(!Save("@X: RGB Image",src)) cerr<<"Could not show the RGB file"<<endl;	
 		
 		//MOUTH REGION DETECTION STEP
 		//Perform Geometric Normalisation and Get Affine Parameters for ROI DETECTION
 		ImageRectangleC facerec = GetFaceCoords(qimg,cascade_name); //Face Recognition Using Viola Jones
		ImageRectangleC res; //to store the result of the ROI in geometrically normalised space
		Tuple2C<ImageC<RealRGBValueC>,Affine2dC> aff; //to store the results of g.normalisation
		try
		{
			Tuple2C<Point2dC,Point2dC> eyes = GetEyeLoc(qimg);
			//AFFINE NORMALISATION OF IMAGE WRT EYES
			aff = AffNormWRTEyes(src, eyes.Data1(),eyes.Data2());
			if(!Save("@X:Affine Normalised RGB",aff.Data1())) cerr<<"Could not display affine normalised image"<<endl;
			res = ImageRectangleC(100,aff.Data1().BRow(),20,110);
			if(!Save("@X: Test Mouth Region Affine Detected Image",ImageC<RealRGBValueC>(aff.Data1(),res))) cerr<<"Could not save the affine ROI image file"<<endl;	
		}
		catch (OmniN::ColossusExceptionNoFaceFoundC er)
	  	{
	   	//If eye estimation doesn't work, just return the bottom half of the face clipped by a window
	   	ImageRectangleC bothalf((facerec.TRow() + facerec.BRow())/2,facerec.BRow() - 10,facerec.LCol() + 20, facerec.RCol() - 20);	 	
	   	res = bothalf;
	  	}		
	  	//Compute the inverse projection rectangle for the normal sized image
	  	Vector2dC im_tl = aff.Data2()*Vector2dC(res.TopLeft().Row(),res.TopLeft().Col());
	  	//cout<<"Aff Image TL = "<<res.TopLeft()<<"\t Image TL = "<<im_tl<<endl;
	  	Vector2dC im_br = aff.Data2()*Vector2dC(res.BottomRight().Row(),res.BottomRight().Col());	  	
	  	//cout<<"Aff Image BR = "<<res.BottomRight()<<"\t Image BR = "<<im_br<<endl;
	  	ImageRectangleC full_img_window(im_tl.Row(),im_br.Row(),im_tl.Col(),im_br.Col());
	  	ImageRectangleC botfacehalf(res);
		ImageC<RealRGBValueC> test_im2(src, full_img_window);
		if(!Save("@X: Test Mouth Region Detected Image",ImageC<RealRGBValueC>(src,full_img_window))) 
		cerr<<"Could not save the ROI image file"<<endl;	
	  	//cout<<"Obtained MOUTH REGION ROI"<<endl;
	  	//Now botfacehalf contains the ROI in g.normed space and full_img_window contains the ROI in image space
		DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > colour_trends; //to store res of clustering
		ImageC<RealT> ppl_img;
		//COLOUR CONVERSION AND ROBUST CLUSTERING STEP
		switch(pix_type)
	  	{
  			case 0:
  			{
  				//Get geometrically normalised image first for robust clustering
  				ImageC<RealRGBValueC> m_img(GetRGBVals(aff.Data1()),botfacehalf);
  				//Perform clustering
  				colour_trends = MCDClustering(m_img,clust_type,h);
  				//Get converted image in normal image space
  				ImageC<RealRGBValueC> m_img2(GetRGBVals(src),full_img_window);
				ppl_img = ProbabilisticPixelLabellingStep(m_img2,colour_trends.Nth(0).Data1(),colour_trends.Nth(1).Data1(), lab_type);
  				break;
  			}
 			case 1:
  			{
  				ImageC<VectorC> m_img(GetNormRGVals(aff.Data1()),botfacehalf);
  				//Perform clustering
  				colour_trends = MCDClustering(m_img,clust_type,h);
  				ImageC<VectorC> m_img2(GetNormRGVals(src),full_img_window);
 				ppl_img = ProbabilisticPixelLabellingStep(m_img2,colour_trends.Nth(0).Data1(),colour_trends.Nth(1).Data1(), lab_type); 
 				break;
  			}
  			case 2:
  			{
  				ImageC<RealHSVValueC> m_img_one(GetHSVVals(aff.Data1()),botfacehalf);
  				ImageC<VectorC> m_img = ConvertToVecImage(m_img_one);
  				//Perform clustering
  				colour_trends = MCDClustering(m_img,clust_type,h);
  				ImageC<RealHSVValueC> m_img2_one(GetHSVVals(src),full_img_window);
  				ImageC<VectorC> m_img2 = ConvertToVecImage(m_img2_one);
 				ppl_img = ProbabilisticPixelLabellingStep(m_img2,colour_trends.Nth(0).Data1(),colour_trends.Nth(1).Data1(), lab_type);				
 				break;
  			}
  			case 3:
  			{
  				ImageC<RealYUVValueC> m_img_one(GetYUVVals(aff.Data1()),botfacehalf);
 				ImageC<VectorC> m_img = ConvertToVecImage(m_img_one);
  				//Perform clustering
  				colour_trends = MCDClustering(m_img,clust_type,h);
  				ImageC<RealYUVValueC> m_img2_one(GetYUVVals(src),full_img_window);
  				ImageC<VectorC> m_img2 = ConvertToVecImage(m_img2_one);
  				ppl_img = ProbabilisticPixelLabellingStep(m_img2,colour_trends.Nth(0).Data1(),colour_trends.Nth(1).Data1(), lab_type);				
	   			break;
  			}
 			case 4:
  			{
  				ImageC<VectorC> m_img(GetCIELabVals(aff.Data1()),botfacehalf);
  				//Perform clustering
  				colour_trends = MCDClustering(m_img,clust_type,h);
  				ImageC<VectorC> m_img2(GetCIELabVals(src),full_img_window);
  				ppl_img = ProbabilisticPixelLabellingStep(m_img2,colour_trends.Nth(0).Data1(),colour_trends.Nth(1).Data1(), lab_type);				
				break;
			}
			case 5:
			{
				ImageC<VectorC> m_img(GetPseudoHueVals(aff.Data1()),botfacehalf);
  				//Perform clustering
  				colour_trends = MCDClustering(m_img,clust_type,h);
  				ImageC<VectorC> m_img2(GetPseudoHueVals(src),full_img_window);
  				ppl_img = ProbabilisticPixelLabellingStep(m_img2,colour_trends.Nth(0).Data1(),colour_trends.Nth(1).Data1(), lab_type);				
				break;
			}
			default:
			{
				ImageC<RealRGBValueC> m_img(GetRGBVals(aff.Data1()),botfacehalf);
  				//Perform clustering
  				colour_trends = MCDClustering(m_img,clust_type,h);
  				ImageC<RealRGBValueC> m_img2(GetRGBVals(src),full_img_window);
  				ppl_img = ProbabilisticPixelLabellingStep(m_img2,colour_trends.Nth(0).Data1(),colour_trends.Nth(1).Data1(), lab_type);			
  				break;		
			}
		}	
		if(!Save("@X: Probabilistic Labelled Mouth Region Image",ppl_img)) 
			cerr<<"Could not save the ROI image file"<<endl;	
		ImageC<UIntT> lbl_img = LabelImage(ppl_img); //This should return a black and white image with lip marked as white and skin as zero
		if(!Save("@X: Segmentation Map of Mouth Region Image",lbl_img)) 
			cerr<<"Could not save the ROI image file"<<endl;	
		//Perform Spatial Cluster Grouping
		//For the moment just with Connected components Analysis
		//ImageC<UIntT> res_im = PerformSpatialRegionIdentification(lbl_img, src.Frame(), reg_type);
		//if(!Save("@X: Binary Lip Image",res_im)) 
		//	cerr<<"Could not save the ROI image file"<<endl;	
		
		ImageC<RealRGBValueC> lip_img(src); ImageC<UIntT> lip_bin(src.Frame(),0);
		RealRGBValueC lip_color(255,10,10);
		for(Array2dIterC<UIntT> it(lbl_img); it; it++)
		{
			if((*it) == 255)
			{
				Index2dC it_ind = it.Index();
				lip_img[it_ind] = lip_color;
				lip_bin[it_ind] = 255;
			}
		}
		if(!Save("@X: Binary Lip Image",lip_img)) 
			cerr<<"Could not save the ROI image file"<<endl;	
		
		if(!Save("@X: Segmented Lip Image",lip_img)) 
			cerr<<"Could not save the ROI image file"<<endl;	
		
		//SCORING CODE GOES IN HERE
		//CURRENT METHOD OF SCORING IS TO COMPARE TWO GROUNT TRUTH IMAGES AND COMPUTE IMAGE SEGMENTATION QUALITY :)
  		//Load up Ground Truth Image
  		ImageC<UIntT> gt_orig;
  		if(!Load(gtimg, gt_orig)) cerr<<"Could not load original Ground Truth Image"<<endl;
  		ImageC<UIntT> gt(gt_orig,full_img_window);
  		if(!Save("@X:GT Image",gt)) cerr<<"Could not save Ground Truth Image"<<endl;
  		RealT score = ComputeSegmentationQuality(lbl_img, gt);
  		cout<<score<<endl;
  		scorelist.Append(score);
		//Now perform lip boundary search to see the spline overlaid
/*
		RealRGBValueC color(10,200,10); RealRGBValueC color2(10,10,200);
		UIntT search = 10;
		ImageC<RealRGBValueC> first_src(src.Copy());
		Array1dC<Point2dC> lip_pt = GetLipContourPoints(lind, 14, lbl_img, search, ppl_img);  	
		for(Array1dIterC<Point2dC> it(lip_pt); it; it++)
	  	{
  			if(first_src.Frame().Contains((*it)))
  				DrawCross(first_src,color,(*it),4);//Draw Contour Points
  		}
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
*/
/*	
  		//Save results	
  		//if(!Save(oimg_seg,lip_img)) cerr<<"Could not save the lip segmentation image"<<endl;
  		//if(!Save(oimg_spl,first_src)) cerr<<"Could not save the lip segmentation image"<<endl;
  		
		//Sort the area
		//SegmentationC seg_map(lbl_img,2);
		//if(!Save("@X:Seg Analysis", seg_map.RandomTaintImage())); 		
*/
		//cin.get();
 	}
	cout<<scorelist<<endl;
	return 0;
}

//Method for fast output
template <typename CL>
DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > MCDClustering(const ImageC<CL> &im, const UIntT &c_type, const RealT &hval)
{
	DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > result;
	switch(c_type)
	{
		case 0:
		{
			//RealT h = 0.0;
			//cout<<"You have chosen to use the Cascade MCD method, please enter the value of h you would like"<<endl;
			//cin>>h;
			//CascadedMCD mcd(LoadImage(im, im_type,irec),h);//load with image type we would like
			CascadedMCD mcd(im,hval);
			//Get the data we want to cluster and input into the clustering algorithm
			result = mcd.Apply();
			break;
		}
		case 1:
		{
			//RealT h = 0.75;
			KSmirnovMCD mcd(im,hval); //load it up with the image type we would like
			//Get the data we want to cluster and input into the clustering algorithm
			result = mcd.Apply();
			break;	
		}
		case 2:
		{	
			UIntT ncl = 3; //declare number of clusters as being equal to 3
			DesignKMeansC km(ncl);
			SampleC<VectorC> smpl = GetSample(im);
			ClassifierC cl = km.Apply(smpl);
			ImageC<UIntT> clustimg(im.Frame(),0);
			UIntT dim = im[0][0].Size();
			DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > out;
			SArray1dC<DListC<Tuple2C<VectorC,Index2dC> > > pix(ncl);
			for(Array2dIterC<CL> it(im); it; it++)
			{
				VectorC v(dim);
				for(UIntT i = 0; i < dim; i++)
				{
					v[i] = (*it)[i];
				}
				UIntT label = cl.Classify(v);
				clustimg[it.Index()] = label;
				//Insert VectorC and Index into the output cluster
				Tuple2C<VectorC,Index2dC> tpl(v.Copy(), it.Index());
				pix[label].Append(tpl);
			}
			//Now obtain the MeanCovariance Values for each cluster
			for(SArray1dIterC<DListC<Tuple2C<VectorC,Index2dC> > > it(pix); it; it++)
			{
				SumsNd2C sum(dim);
				bool SampleStats = false;
				for(DLIterC<Tuple2C<VectorC,Index2dC> > it2((*it)); it2; it2++)
				{
					sum += (*it2).Data1();
				}
				Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > elem(sum.MeanCovariance(SampleStats).Copy(),(*it));
				out.Append(elem.Copy());
			}
			result = out.Copy();			
		}
		case 3:
		{
			UIntT ncl = 3; //declare number of clusters as being equal to 3
			DesignFuzzyCMeansClusterC km(ncl);
			SampleC<VectorC> smpl = GetSample(im);
			ClassifierC cl = km.Apply(smpl);
			ImageC<UIntT> clustimg(im.Frame(),0);
			UIntT dim = im[0][0].Size();
			DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > out;
			SArray1dC<DListC<Tuple2C<VectorC,Index2dC> > > pix(ncl);
			for(Array2dIterC<CL> it(im); it; it++)
			{
				VectorC v(dim);
				for(UIntT i = 0; i < dim; i++)
				{
					v[i] = (*it)[i];
				}
				UIntT label = cl.Classify(v);
				clustimg[it.Index()] = label;
				//Insert VectorC and Index into the output cluster
				Tuple2C<VectorC,Index2dC> tpl(v.Copy(), it.Index());
				pix[label].Append(tpl);
			}
			//Now obtain the MeanCovariance Values for each cluster
			for(SArray1dIterC<DListC<Tuple2C<VectorC,Index2dC> > > it(pix); it; it++)
			{
				SumsNd2C sum(dim);
				bool SampleStats = false;
				for(DLIterC<Tuple2C<VectorC,Index2dC> > it2((*it)); it2; it2++)
				{
					sum += (*it2).Data1();
				}
				Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > elem(sum.MeanCovariance(SampleStats).Copy(),(*it));
				out.Append(elem.Copy());
			}
			result = out.Copy();					
		}
		
		default:
		{
			CascadedMCD mcd(im,hval);//load with image type we would like
			//Get the data we want to cluster and input into the clustering algorithm
			result = mcd.Apply();
			break;
		
		}
	}
	return result;
}

Array1dC<Point2dC> GetLipContourPoints(const DListC<Index2dC> &lp, const UIntT &res, const ImageC<UIntT> &img, const RealT &search, const ImageC<RealT> &jim) 
{
	//cout<<"Binary Image Frame = "<<img.Frame()<<endl;
	//Convert to an SArray1dC<Index2dC> for ellipse fitting
	SArray1dC<Point2dC> samp(lp.Size());
	UIntT ind = 0;
	for(DLIterC<Index2dC> d(lp); d; d++)
	{
		Point2dC pt((*d).Row(),(*d).Col());
		samp[ind++] = pt.Copy();
	}
	Ellipse2dC ell;
	FitEllipse(samp,ell);
	//Now Generate Points on the Ellipse
	RealT inc = (2*pi)/(RealT)res;
	Array1dC<Point2dC> out(res);
	for(UIntT i = 0; i < res; i++)
	{
		RealT angle = (RealT)i * inc;
		out[i] = ell.Point(angle).Copy();
	}
	//cout<<"Ellipse Point Generation"<<endl;
	Array1dC<Point2dC> first(1); first[0] = out[0];
	out.Append(first);
	//Using generated points, Get B-Spline
	UIntT order = 3; UIntT ncp = 11;
	BSplineC bspl(order,ncp, BSplineC::UOPEN);
	Array1dC<Point2dC> cpts = bspl.CalculateControlPoints(out, BSplineC::CHORDLENGTH, BSplineC::UOPEN,order, ncp);	
	//cout<<"Calculated Control Points"<<cpts<<endl;
	Array1dC<Point2dC> lip_ctr(res); UIntT ctr_ind = 0;
	RealT incr = 1.0/(RealT)res;
	for(RealT pval = 0; pval < res; pval ++)		
	{
		RealT par = (RealT)pval*incr;
		LinePP2dC normal = bspl.GetNormalLine(par); //Get B-Spline Normal Line
		Point2dC start = normal.Point(0.0);
		lip_ctr[ctr_ind] = start.Copy();
		//Now just search from the start to the end for the binary condition to end
		for(RealT i = -5; i < search; i++)
		{
			Point2dC cur = normal.Point(i); Point2dC next = normal.Point(i+1);
			if(jim.Frame().Contains(next))
			{
				if((jim[cur] >= 0)&&(jim[next] < 0))
				{
					lip_ctr[ctr_ind] = normal.Point(i).Copy();
				}
			}
		}
		ctr_ind++;
	}
	//cout<<"Exited Loop "<<lip_ctr<<endl;
	return lip_ctr;
}

//CentralColumnRegions()
//This method assigns all the regions containing a non-central column to be part of the background
SegmentationC CentralColumnRegions(const SegmentationC &seg_in)
{
	SegmentationC seg = seg_in;
	//Get the central column of the image inside the segmentation object
	IndexC cent_col = seg.SegMap().Frame().Center().Col();
	//if object has only one label, then can't do anything to it and return
	if(seg.Labels() <= 1)
	      return seg;
	SArray1dC<Tuple2C<IndexRange2dC,UIntT> > bndarea = seg.BoundsAndArea();
    SArray1dC<UIntT> area = seg.Areas();
	// Assign new labels to the regions according to whether they contain the central image column or not
    IntT newLabel = 1;
	
	SArray1dIterC<Tuple2C<IndexRange2dC,UIntT> > it1(bndarea);
	SArray1dIterC<UIntT> it(area);
    *it = 0;
    it.Next();
	it1.Next();
    for(;it;it++) 
    {
		if (!(*it1).Data1().Contains(Index2dC((*it1).Data1().Center().Row(),cent_col))) 
			*it = 0;
		else 
			*it = newLabel++;
		
		it1++;
    }
    // Remove small components
	
    for(Array2dIterC<UIntT> iti(seg.SegMap());iti;iti++)
		*iti = area[*iti];
    SegmentationC res(seg.SegMap(),newLabel);
    return res;    
}

//ProbabilisticPixelLabellingStep()
// This performs probabilistic pixel labelling to produce a binary image based on the cost function specified (2 are available)
template <typename T>
ImageC<RealT> ProbabilisticPixelLabellingStep(const ImageC<T> &im, const MeanCovarianceC &sk, const MeanCovarianceC &lp, const UIntT &type)
{
	ImageC<RealT> res;
	switch(type)
	{
	case 0: //This uses the normal J function
	{	
		res = ProbabilisticPixelLabelling(im,sk,lp);	
		break;
	}
	case 1:
	{
		res = ProbabilisticPixelLabelling(im,sk,MeanCovarianceC(sk.Number(),lp.Mean(),sk.Covariance()));	
		break;	
	}
	default:
	{
		res = ProbabilisticPixelLabelling(im,sk,lp);	
		break;
	}
	}
	return res;
}

//GetLargestConnectedComponent()
//Obtain largest, non-zero connected component given a binary image
ImageC<UIntT> GetLargestConnectedComponent(const ImageC<UIntT> &im, const ImageRectangleC &imrec)
{
	
	ConnectedComponentsC<UIntT> connected; 
	Tuple2C<ImageC<UIntT>, UIntT> result = connected.Apply(im) ; 
	if(!Save("@X:CC Analysis", SegmentationC(result).RandomTaintImage())) 
		cerr << "\n Error failed to load output image " ;	
	//SegmentationC sm(result);
	SegmentationC sm = CentralColumnRegions(SegmentationC(result));
	if(!Save("@X: Central Column Stuff",sm.RandomTaintImage())) cerr<<"Could not save central column image"<<endl;
	if(sm.Labels() <= 1 )
	return sm.SegMap();
	
	SArray1dC<UIntT> sarea = sm.Areas().Copy();
	sarea.Sort();
	IntT lar = sarea[1];
	UIntT label = 0;
	for(IntT i = 0; i < (IntT)sm.Areas().Size(); i++)
	{
		if((sm.Areas())[i] == (UIntT)lar)
			label = i;
	}
	//In the original image, label the lip pixels as well
	//if(!Save("@X:JCC", lip_img)); 
	DListC<Index2dC> lind;
	//ImageC<RealRGBValueC> lip_img(src.Copy());
	ImageC<UIntT> res(imrec, 0);
	//RealRGBValueC lip_label(255,50,50);
	UIntT lip_label = 255;
	for(Array2dIterC<UIntT> it(result.Data1()); it; it++)
	{
		if((*it)== label)
		{
			lind.Append(it.Index().Copy());
			res[it.Index()] = lip_label;
		}
	}
	return res;
}

//PerformSpatialRegionIdentification()
//This performs spatial-property based region identification depending on the type specified
ImageC<UIntT> PerformSpatialRegionIdentification(const ImageC<UIntT> &im, const ImageRectangleC &imrec, const UIntT &type)
{
		ImageC<UIntT> res;
		switch(type)
		{
			case 0:
			{
				res = GetLargestConnectedComponent(im, imrec);
				break;
			}
			default:
			{
				PerformSpatialRegionIdentification(im,imrec,0); //default is recursive set to option 0 to preserve code consistency
			}
		}
		return res;
}
