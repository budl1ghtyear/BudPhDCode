#include "LipSegmentation.hh"

LipSegmentationC::LipSegmentationC(const ImageC<RealRGBValueC> &img)
{
		inp_img = img.Copy();
		inp_img_name = "UNKNOWN";
}

LipSegmentationC::LipSegmentationC(const FilenameC &file)
{
	inp_img_name = file;
	if(!Load(file,inp_img)) cerr<<"Could not load the file, INVALID FILE NAME? "<<endl;
}

ColourSpaceImage LipSegmentationC::ColourSpaceConversion(const PixelType &ptype)
{
	ImageC<VectorC> m_img;
	switch(ptype)
	{
		case RGB:
		{
			m_img = ConvertToVecImage(inp_img);
			break;
		}
		case NormRG:
		{
			m_img = GetNormRGVals(inp_img);
			break;
		}
		case HSV:
		{
 			ImageC<RealHSVValueC> m_img_one(GetHSVVals(inp_img));
  			m_img = ConvertToVecImage(m_img_one);
			break;
		}
		case YUV:
		{
 			ImageC<RealYUVValueC> m_img_one(GetYUVVals(inp_img));
 			m_img = ConvertToVecImage(m_img_one);
			break;
		}
		case CIE:
		{  				
			m_img = GetCIELabVals(inp_img);	
			break;				
		}
		case PHue:
		{
			m_img = GetPseudoHueVals(inp_img);
			break;
		}
		default:
		{
			m_img = ConvertToVecImage(inp_img);
			break;
			//default is RGB
		}
	}
	ColourSpaceImage res(m_img,ptype);
	return res;
}

Tuple4C<Affine2dC, ImageRectangleC, ImageRectangleC,ImageC<VectorC> > LipSegmentationC::AffineMouthRegionDetection(const ColourSpaceImage &im)
{
	FilenameC cascade_name = "/opt/share/opencv/haarcascades/haarcascade_frontalface_alt.xml"; //REMEMBER TO CHANGE THIS ON DIFFERENT INSTALLATIONS
	//FilenameC cascade_name = "/usr/share/opencv/haarcascades/haarcascade_frontalface_alt.xml"; //REMEMBER TO CHANGE THIS ON DIFFERENT INSTALLATIONS
	ImageRectangleC facerec = GetFaceCoords(inp_img_name,cascade_name); //perform face detection using OpenCV Viola-Jones
	ImageC<RealRGBValueC> face_img(inp_img,facerec);
	//if(!Save("@X:Face Image",face_img)) cerr<<"Face cannot be shown"<<endl;
	ImageRectangleC res; //to store result
	Tuple2C<ImageC<VectorC>,Affine2dC> aff; //to store the results of g.normalisation
	try
	{
		Tuple2C<Point2dC,Point2dC> eyes = GetEyeLoc(inp_img); //Method is in MouthRegion.hh
		//cout<<"eyes "<<eyes<<endl;
		//AFFINE NORMALISATION OF IMAGE WRT EYES
		aff = AffNormWRTEyes(im.GetImage(), eyes.Data1(),eyes.Data2()); //Method is in MouthRegion.hh
		res = ImageRectangleC(100,aff.Data1().BRow(),20,110);
		#ifdef DEBUGSHARED
		//if(!Save("@X:Affine Normalised RGB",aff.Data1())) cerr<<"Could not display affine normalised image"<<endl;
		if(!Save("@X: Test Mouth Region Affine Detected Image",ImageC<RealRGBValueC>(inp_img,res))) cerr<<"Could not save the affine ROI image file"<<endl;	
		#endif
	}
	catch (OmniN::ColossusExceptionNoFaceFoundC er)
	{
	   	//If eye estimation doesn't work, just return the bottom half of the face clipped by a window
	   	ImageRectangleC bothalf((facerec.TRow() + facerec.BRow())/2,facerec.BRow() - 10,facerec.LCol() + 20, facerec.RCol() - 20);	 	
	   	res = bothalf;
		//Instantiate UNIT AFFINE transform (i.e. one that just transforms the image onto itself)
		Affine2dC aff_unit(Vector2dC(1.0,1.0),(RealT)0.0,(RealT)0.0,Vector2dC(0.0,0.0));
		aff.Data1() = im.GetImage().Copy();
		aff.Data2() = aff_unit;		
	 }		
	 //Compute the inverse projection rectangle for the normal sized image
	 Vector2dC im_tl = aff.Data2()*Vector2dC(res.TopLeft().Row(),res.TopLeft().Col());
	 Vector2dC im_br = aff.Data2()*Vector2dC(res.BottomRight().Row(),res.BottomRight().Col());	  	
	 ImageRectangleC full_img_window(im_tl.Row(),im_br.Row(),im_tl.Col(),im_br.Col());
	 ImageRectangleC botfacehalf(res);		
	 #ifdef DEBUGSHARED
 	 ImageC<RealRGBValueC> test_im2(inp_img, full_img_window);
	 if(!Save("@X: Test Mouth Region Detected Image",ImageC<RealRGBValueC>(inp_img,full_img_window))) cerr<<"Could not save the ROI image file"<<endl;	
	 cout<<"Aff Image BR = "<<res.BottomRight()<<"\t Image BR = "<<im_br<<endl;
	 cout<<"Aff Image TL = "<<res.TopLeft()<<"\t Image TL = "<<im_tl<<endl;	 
	 cout<<"Obtained MOUTH REGION ROI"<<endl;
	 #endif
	 aff_transform = aff.Data2();
	 //Now botfacehalf contains the ROI in g.normed space and full_img_window contains the ROI in image space
	 return (Tuple4C<Affine2dC, ImageRectangleC, ImageRectangleC,ImageC<VectorC> >(aff.Data2(),botfacehalf,full_img_window,ImageC<VectorC>(aff.Data1(),botfacehalf)));
}

Tuple2C<MeanCovarianceC,MeanCovarianceC> LipSegmentationC::LipClusteringStep(const ClusteringType &type, const ImageC<VectorC> &im, const RealT &hval=(RealT)0.75, const UIntT &nclust=(UIntT)3)
{
		//if(!Save("@X:CLI Img",im.GetImage())) cerr<<"CLI Image Not Shown"<<endl;
		LipClusteringC<VectorC> clust(im,hval,nclust);
		DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > trends = clust.Apply(type);	//Because of the way the clustering works, this does not need to be sorted by size
		color_trends = Tuple2C<MeanCovarianceC,MeanCovarianceC>(trends.PopFirst().Data1(), trends.PopFirst().Data1());
		//Draw clustered image
		//cout<<"Num colour trends extracted - "<<trends.Size()<<endl;
		ImageC<UIntT> clust_img(im.Frame(),0); UIntT label = 0;
		for(DLIterC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > it(trends); it; it++)
		{
			for(DLIterC<Tuple2C<VectorC,Index2dC> > it2((*it).Data2()); it2; it2++)
			{
				clust_img[(*it2).Data2()] = label;
			}
			label++;
		}
		SegmentationC seg(clust_img,label); 
		return (color_trends);
}

ImageC<UIntT> LipSegmentationC::ProbabilisticPixelLabellingStep(const ImageC<VectorC> &im, const MeanCovarianceC &sk, const MeanCovarianceC &lp, const LabellingType &type)
{
	ImageC<RealT> res;
	switch(type)
	{
		case JNormal: //This uses the normal J function
		{	
			res = ProbabilisticPixelLabelling(im,sk,lp);//Method is in CostFunctions.hh
			break;
		}
		case JSkin:
		{
			res = ProbabilisticPixelLabelling(im,sk,MeanCovarianceC(sk.Number(),lp.Mean(),sk.Covariance())); //Method is in CostFunctions.hh	
			break;	
		}
		default:
		{
			res = ProbabilisticPixelLabelling(im,sk,lp);//Method is in CostFunctions.hh	
			break;
		}
	}
	//if(!Save("@X:PPL Image",res)) cerr<<"Could not save the labelled image"<<endl;
	//if(!Save("@X:PPL Binary Image",LabelImage(res))) cerr<<"Could not save the labelled image"<<endl;
	return (LabelImage(res)); //Method LabelImage() is in CostFunctions.hh
}

ImageC<UIntT> LipSegmentationC::RegionIdentificationStep(const ImageC<UIntT> &im, const RegionIdentificationType &type)
{
	ImageC<UIntT> res;
	switch(type)
	{
		case CC: //for largest connected component
		{
			ConnectedComponentsC<UIntT> connected; 
			Tuple2C<ImageC<UIntT>, UIntT> result = connected.Apply(im) ; 
			#ifdef DEBUGSHAREDSHARED
				if(!Save("@X:CC Analysis", SegmentationC(result).RandomTaintImage())) cerr << "\n Error failed to load output image " ;	
			#endif
			SegmentationC sm(result);
			if(sm.Labels() <= 1 )
			{
				res = sm.SegMap().Copy();			
				break;
			}
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
			ImageC<UIntT> imres(im.Frame(), 0);
			//RealRGBValueC lip_label(255,50,50);
			UIntT lip_label = 255;
			for(Array2dIterC<UIntT> it(result.Data1()); it; it++)
			{
				if((*it)== label)
				{
					lind.Append(it.Index().Copy());
					imres[it.Index()] = lip_label;
				}
			}
			res = imres.Copy();
			break;
		}
		case CCentral: //for largest connected component along central image axis
		{
			ConnectedComponentsC<UIntT> connected; 
			Tuple2C<ImageC<UIntT>, UIntT> result = connected.Apply(im) ; 
			#ifdef DEBUGSHARED
				if(!Save("@X:CC Analysis", SegmentationC(result).RandomTaintImage())) cerr << "\n Error failed to load output image " ;	
			#endif
			SegmentationC sm = CentralColumnRegions(SegmentationC(result)); //Method is in CostFunctions.hh
			if(sm.Labels() <= 1 )
			{
				res = sm.SegMap().Copy();			
				break;
			}
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
			ImageC<UIntT> imres(im.Frame(), 0);
			//RealRGBValueC lip_label(255,50,50);
			UIntT lip_label = 255;
			for(Array2dIterC<UIntT> it(result.Data1()); it; it++)
			{
				if((*it)== label)
				{
					lind.Append(it.Index().Copy());
					imres[it.Index()] = lip_label;
				}
			}
			res = imres.Copy();
			break;			
		}
		default:
		{
			res = RegionIdentificationStep(im, (RegionIdentificationType)0).Copy(); //default case is with largest connected component
			break;
		}
	}
	return res; //PLEASE NOTE THAT THIS IMAGE IS JUST THE MOUTH_REGION IMAGE AND THE LIP VALUE IS 255
}

Tuple3C<ImageC<UIntT>, MeanCovarianceC, MeanCovarianceC> LipSegmentationC::Apply(const PixelType &ptype,const ClusteringType &ctype,const LabellingType &ltype,const RegionIdentificationType &rtype,const RealT &hval=(RealT)0.75, const UIntT &nclust=(UIntT)3)
{
	//Just perform the modules one by one :) 
	ColourSpaceImage col_img(ColourSpaceConversion(ptype));
	#ifdef DEBUGSHARED
		//Check that we have something
		cout<<"Finished Colour Space Conversion"<<endl;
		cout<<"Colour Space is - "<<col_img.GetImageType()<<" and the image is of size - "<<col_img.GetImage().Frame()<<endl;		
	#endif
	Tuple4C<Affine2dC, ImageRectangleC, ImageRectangleC,ImageC<VectorC> > amri = AffineMouthRegionDetection(col_img);
	#ifdef DEBUGSHARED
		cout<<"Finished Affine Geometric Normalisation and Mouth Region Identification"<<endl;
		cout<<"Affine Transform is "<<amri.Data1()<<" with Affine Mouth Frame as "<<amri.Data2()<<" and original Image Mouth Frame is "<<amri.Data3()<<endl;
	#endif
	Tuple2C<MeanCovarianceC,MeanCovarianceC> col_trends = LipClusteringStep(ctype,ImageC<VectorC>(amri.Data4(),amri.Data2()),hval,nclust);
	#ifdef DEBUGSHARED
		cout<<"Finished Clustering"<<endl;
		cout<<"Skin estimate - "<<col_trends.Data1()<<" Lip estimate = "<<col_trends.Data2()<<endl;
	#endif	
	ImageC<UIntT> ppl_img = ProbabilisticPixelLabellingStep(ImageC<VectorC>(col_img.GetImage(),amri.Data3()) , col_trends.Data1(), col_trends.Data2(), ltype);
	#ifdef DEBUGSHARED
		cout<<"Finished Pixel Labelling"<<endl;
	#endif
	ImageC<UIntT> rgn_img = RegionIdentificationStep(ppl_img,rtype);
	#ifdef DEBUGSHARED
		cout<<"Finished Region Identification"<<endl;
	#endif
	//NOW WE HAVE A BINARY LIP IMAGE, WE CAN RETURN IT
	return (Tuple3C<ImageC<UIntT>, MeanCovarianceC, MeanCovarianceC>(rgn_img,col_trends.Data1(), col_trends.Data2()));	
}
