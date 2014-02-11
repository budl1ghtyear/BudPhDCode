#include "SegmentationC.hh"

SegmentationC::SegmentationC(const FilenameC &file, PixelType &ptype, ClusteringType &ctype, LabellingType &ltype, RegionIdenticationType &rtype, const RealT &hvalue)
{
	in_file = file;
	ImageC<RealRGBValueC> src;   
	if(!Load(file, src)) cerr<<"Loading RAVL Image Failed"<<endl;
	#ifdef DEBUG
	if(!Save("@X: RGB Image",src)) cerr<<"Could not show the RGB file"<<endl;		
	#endif
}
//: Constructor defining the Input filename as well as processing options
//: Apply Method Will be called inside the constructor in addition to being defined explicitly for consistency purposes

SegmentationC::SegmentationC(const FilenameC &file, PixelType &ptype, ClusteringType &ctype, LabellingType &ltype, RegionIdenticationType &rtype)
{
	
}
//:Constructor method for default MCD options with default "h" value of 0.75


SegmentationC::SegmentationC(const FilenameC &file, PixelType &ptype, ClusteringType &ctype, LabellingType &ltype, RegionIdenticationType &rtype, const UIntT &num_cl)
{
	
}
//:Constructor with method for instantiating K-Means and Fuzzy C Means algorithms with a predefined number of clusters

Tuple3C<Affine2dC, ImageRectangleC, ImageRectangleC> SegmentationC::AffineMouthRegionDetection(const ImageC<RealRGBValueC> &im, const PixelType &ptype)
{
 	ImageRectangleC facerec = GetFaceCoords(in_file,cascade_name); //Face Recognition Using Viola Jones
	ImageRectangleC res; //to store the result of the ROI in geometrically normalised space
	Tuple2C<ImageC<RealRGBValueC>,Affine2dC> aff; //to store the results of g.normalisation
	try
	{
		Tuple2C<Point2dC,Point2dC> eyes = GetEyeLoc(in_file);
		//AFFINE NORMALISATION OF IMAGE WRT EYES
		aff = AffNormWRTEyes(im, eyes.Data1(),eyes.Data2());
		#ifdef DEBUG
		if(!Save("@X:Affine Normalised RGB",aff.Data1())) cerr<<"Could not display affine normalised image"<<endl;
		#endif
		res = ImageRectangleC(100,aff.Data1().BRow(),20,110);
		#ifdef DEBUG
		if(!Save("@X: Test Mouth Region Affine Detected Image",ImageC<RealRGBValueC>(aff.Data1(),res))) cerr<<"Could not save the affine ROI image file"<<endl;	
		#endif
	}
	catch (OmniN::ColossusExceptionNoFaceFoundC er)
	{
	 	//If eye estimation doesn't work, just return the bottom half of the face clipped by a window
	   	ImageRectangleC bothalf((facerec.TRow() + facerec.BRow())/2,facerec.BRow() - 10,facerec.LCol() + 20, facerec.RCol() - 20);	 	
	   	res = bothalf;
	}		
	//Compute the inverse projection rectangle for the normal sized image
	Vector2dC im_tl = aff.Data2()*Vector2dC(res.TopLeft().Row(),res.TopLeft().Col());
	Vector2dC im_br = aff.Data2()*Vector2dC(res.BottomRight().Row(),res.BottomRight().Col());	  	
	ImageRectangleC full_img_window(im_tl.Row(),im_br.Row(),im_tl.Col(),im_br.Col());
	ImageRectangleC botfacehalf(res);
	ImageC<RealRGBValueC> test_im2(src, full_img_window);
	#ifdef DEBUG
	cout<<"Aff Image TL = "<<res.TopLeft()<<"\t Image TL = "<<im_tl<<endl;
	cout<<"Aff Image BR = "<<res.BottomRight()<<"\t Image BR = "<<im_br<<endl;
	if(!Save("@X: Test Mouth Region Detected Image",ImageC<RealRGBValueC>(src,full_img_window))) cerr<<"Could not save the ROI image file"<<endl;	
	cout<<"Obtained MOUTH REGION ROI"<<endl;
	#endif
	//Now botfacehalf contains the ROI in g.normed space and full_img_window contains the ROI in image space
	Tuple3C<Affine2dC, ImageRectangleC, ImageRectangleC> out(aff.Copy(), botfacehalf, full_img_window);
	return out;	
}
//:Get Affine Normalised Image For Use In Clustering as well as resulting ImageRectangle Windows

Tuple2C<MeanCovarianceC, MeanCovarianceC> SegmentationC::GetChromaticClusters(ColourSpaceImage &im, const ClusteringType &ctype, const RealT &hval)
{
}
//:Get Cluster Trends Info For Probabilistic Pixel Labelling
ImageC<UIntT> SegmentationC:: ProbabilisticPixelLabellingStep(const ImageC<T> &im, const MeanCovarianceC &sk, const MeanCovarianceC &lp, const UIntT &type)
{
}
//:Perform Probabilistic Pixel Labelling
ImageC<UIntT> SegmentationC::RegionIdenticationStep(const ImageC<UIntT> &im, const ImageRectangleC &imrec, const RegionIdenticationType &rtype)
{
}
//:Perform Region Identification
