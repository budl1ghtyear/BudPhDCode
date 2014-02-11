//File: SegmentationDemo.cc
// Author: Bud Goswami
// Input options: Input image directory, pixel reppresentation type, clustering type and region identification type

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

DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > GetColourTrends(const ImageC<RGBValueC> &im, const UIntT &p_type, const UIntT &c_type, const RealT &hval);
template <typename CL>
DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > MCDClustering(const ImageC<CL> &im, const UIntT &c_type, const RealT &hval);
Array1dC<Point2dC> GetLipContourPoints(const DListC<Index2dC> &lp, const UIntT &res, const ImageC<UIntT> &img, const RealT &search, const ImageC<RealT> &jim);
SegmentationC CentralColumnRegions(SegmentationC &seg);

int main(int nargs,char **argv) 
{
	//OBTAIN INPUT OPTIONS
	OptionC opt(nargs,argv);
	DirectoryC qimg_dir = opt.String("i","TestSet/","Input Query Image Directory with ppm images");
	DirectoryC outimg_dir = opt.String("o","Output/","Output Image Directory");
	FilenameC cascade_name = opt.String("c","/usr/share/opencv/haarcascades/haarcascade_frontalface_alt.xml","Path for the Haar Classifier");
	StringC extension = opt.String("e",".ppm","Input Image Extension");
	UIntT pix_type = opt.Int("p",0,"Pixel Output Type - 0=RealRGBValueC, 1=NormalisedRG, 2=RealHSVValueC, 3=RealYUVValueC, 4=CIELab, 5 = PseudoHue");
	UIntT clust_type = opt.Int("l",0,"Robust Clustering Type - 0=CascadedMCD, 1=MCD with KS, 2=K Means, 3= Fuzzy C Means");
	RealT h = opt.Real("h",0.75,"Confidence Value for Methods");
	opt.Check();
	//
	//RUN LOOP FOR DIRECTORY IMAGES
	DListC<StringC> qfile = qimg_dir.FiltList("*"+extension);
	for(DLIterC<StringC> dt(qfile); dt; dt++)
	{
		//GET NECESSARY FILES
		FilenameC qimg = qimg_dir+(*dt);		
		FilenameC oimg_seg = outimg_dir+"SegmentationResults/"+qimg.BaseNameComponent()+"SEG.ppm";
		FilenameC oimg_spl = outimg_dir+"SplineResults/"+qimg.BaseNameComponent()+"SPL.ppm";
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
	  	cout<<"Aff Image TL = "<<res.TopLeft()<<"\t Image TL = "<<im_tl<<endl;
	  	Vector2dC im_br = aff.Data2()*Vector2dC(res.BottomRight().Row(),res.BottomRight().Col());	  	
	  	cout<<"Aff Image BR = "<<res.BottomRight()<<"\t Image BR = "<<im_br<<endl;
	  	ImageRectangleC full_img_window(im_tl.Row(),im_br.Row(),im_tl.Col(),im_br.Col());
	  	ImageRectangleC botfacehalf(res);
		ImageC<RealRGBValueC> test_im2(src, full_img_window);
		if(!Save("@X: Test Mouth Region Detected Image",ImageC<RealRGBValueC>(src,full_img_window))) 
		cerr<<"Could not save the ROI image file"<<endl;	
	  	cout<<"Obtained MOUTH REGION ROI"<<endl;
	  	//Now botfacehalf contains the ROI in g.normed space and full_img_window contains the ROI in image space
		DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > colour_trends; //to store res of clustering
		ImageC<RealT> ppl_img;
		
		//COLOUR CONVERSION AND ROBUST CLUSTERING STEP
		if(!Save("@X: Probabilistic Labelled Mouth Region Image",ppl_img)) 
		cerr<<"Could not save the ROI image file"<<endl;	
		ImageC<UIntT> lbl_img = LabelImage(ppl_img); //This should return a black and white image with lip marked as white and skin as zero
		if(!Save("@X: Segmentation Map of Mouth Region Image",lbl_img)) 
		cerr<<"Could not save the ROI image file"<<endl;	
		//Perform Spatial Cluster Grouping
		//For the moment just with Connected components Analysis
		ConnectedComponentsC<UIntT> connected; 
		Tuple2C<ImageC<UIntT>, UIntT> result = connected.Apply(lbl_img) ; 
		if(!Save("@X:CC Analysis", SegmentationC(result).RandomTaintImage())) 
		//Now remove all the remaining components that Do NOT contain the central column in the image
		//SegmentationC step_one(result);
		


		cerr << "\n Error failed to load output image " ;	
		SegmentationC sm(result);
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
		DListC<Index2dC> lind;
		ImageC<RealRGBValueC> lip_img(src.Copy());
		RealRGBValueC lip_label(255,50,50);
		for(Array2dIterC<UIntT> it(result.Data1()); it; it++)
		{
			if((*it)== label)
			{
				lind.Append(it.Index().Copy());
				lip_img[it.Index()] = lip_label;
			}
		}
		if(!Save("@X:JCC", lip_img)); 
 	}
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
	cout<<"Binary Image Frame = "<<img.Frame()<<endl;
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
/*
SegmentationC CentralColumnRegions(SegmentationC &seg)
{
	if(seg.Labels() <= 1)
	      return seg;
       	SArray1dC<UIntT> area = seg.Areas();
    	UIntT thrSize = 100;
    	// Assign new labels to the regions according their sizes and
    	// another requirements.
    	IntT newLabel = 1;
    	SArray1dIterC<UIntT> it(area);
    	*it = 0;
    	it.Next();
    	for(;it;it++) 
    	{
    	  if (*it < ((UIntT) thrSize)) 
		*it = 0;
    	  else 
		*it = newLabel++;
    	}
    	
    	// Remove small components
    	for(Array2dIterC<UIntT> iti(seg);iti;iti++)
    	  *iti = area[*iti];
    	
    	//seg.Labels() = newLabel;
    	return seg;    
}
*/
//GetColourTrends() Function
//Inputs - Clipped Mouth-Region ROI as an ImageC<RealRGBValueC>, pixel-representation opyion, clustering option, value of h-confidence to use

