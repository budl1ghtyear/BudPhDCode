//CostFunctions.cc
//Author - Bud Goswami
//Date - 18.03.09
//Defines the cost functions and scoring related functions for cluster segmentation and grouping

#include "CostFunctions.hh"


bool IsInsideEllipse(const Ellipse2dC &ell, const Index2dC &ind)
{
 /*  I assume you also know the location of the ellipse's center. Call that (x0,y0).
Let t be the counterclockwise angle the major axis makes with respect to the
x-axis. Let a and b be the semi-major and semi-minor axes, respectively. If
P = (x,y) is an arbitrary point then do this:

 X = (x-x0)*cos(t)+(y-y0)*sin(t); % Translate and rotate coords.
 Y = -(x-x0)*sin(t)+(y-y0)*cos(t); % to align with ellipse

If

 X^2/a^2+Y^2/b^2

is less than 1, the point P lies inside the ellipse. If it equals 1, it is right on
the ellipse. If it is greater than 1, P is outside.*/

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

bool mscomp(const Tuple2C<MeanCovarianceC, DListC<Tuple2C<VectorC,Index2dC> > > &dat1, const Tuple2C<MeanCovarianceC, DListC<Tuple2C<VectorC,Index2dC> > > &dat2)
{
	return (dat1.Data2().Size() >= dat2.Data2().Size() ? true:false);
	//return (dat1.Data2().Size() <= dat2.Data2().Size() ? true:false);
}

DListC<Index2dC> JCC(DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > &clusts,const ImageC<RealRGBValueC> &img)
{
	//oimg needs to be the RGB normalised image
	//Show off the clustered Image
	UIntT lab = 0;
	ImageC<UIntT> m_img(img.Frame(),0);
	for(DLIterC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > it(clusts); it; it++)
	{
		for(DLIterC<Tuple2C<VectorC,Index2dC> > it2((*it).Data2()); it2; it2++)
		{
			m_img[(*it2).Data2()] = lab;
		}
		lab++;
	}
	SegmentationC dummyseg(m_img, lab);
	if(!Save("@X: Test Segmentation Image",dummyseg.RandomImage()))
		cerr<<"Failed to show segmentation image"<<endl;
	
	
	clusts.MergeSort(mscomp);
	ImageC<IntT> im(img.Frame());
	im.Fill(0);
	MeanCovarianceC scl = clusts.Nth(0).Data1();
	MeanCovarianceC lcl = clusts.Nth(1).Data1();
	ImageC<UIntT> dummy(img.Frame());
	dummy.Fill(0);
	RealT sfirst = -Log(scl.Covariance().Det());
	RealT lfirst = -Log(lcl.Covariance().Det());
	for(Array2dIter2C<UIntT, RealRGBValueC> it(dummy, img); it; it++)
	{
		VectorC v(it.Data2().Red(), it.Data2().Green(), it.Data2().Blue());
		RealT Jlips = lfirst - lcl.MahalanobisDistance(v);
		RealT Jskin = sfirst - scl.MahalanobisDistance(v);
		RealT j = Jlips - Jskin;
		if(j > 0)
			it.Data1() = 1;
	}	
	ConnectedComponentsC<UIntT> connected ; 
	// Apply the algorithm 
	Tuple2C<ImageC<UIntT>, UIntT> result = connected.Apply(dummy);
	SegmentationC sm(result);
	SArray1dC<UIntT> sarea = sm.Areas().Copy();
	sarea.Sort();
	//Look for the label corresponding to the second largest area (after skin)
	IntT lar = sarea[1];
	UIntT label = 0;
	for(IntT i = 0; i < (IntT)sm.Areas().Size(); i++)
	{
		if((sm.Areas())[i] == (UIntT)lar)
			label = i;
	}
	//Insert all label pixels to be lips
	DListC<Index2dC> lind;
	for(Array2dIterC<UIntT> it(result.Data1()); it; it++)
	{
		if((*it)== label)
			lind.Append(it.Index().Copy());
	}
	return lind;	
}

DListC<Index2dC> JCC(DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > &clusts,const ImageC<TFVectorC<RealT,3> > &img)
{
	//oimg needs to be the RGB normalised image
	//Show off the clustered Image
	UIntT lab = 0;
	ImageC<UIntT> m_img(img.Frame(),0);
	for(DLIterC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > it(clusts); it; it++)
	{
		for(DLIterC<Tuple2C<VectorC,Index2dC> > it2((*it).Data2()); it2; it2++)
		{
			m_img[(*it2).Data2()] = lab;
		}
		lab++;
	}
	SegmentationC dummyseg(m_img, lab);
	if(!Save("@X: Test Segmentation Image",dummyseg.RandomImage()))
		cerr<<"Failed to show segmentation image"<<endl;
	clusts.MergeSort(mscomp);
	ImageC<IntT> im(img.Frame());
	im.Fill(0);
	MeanCovarianceC scl = clusts.Nth(0).Data1();
	MeanCovarianceC lcl = clusts.Nth(1).Data1();
	ImageC<UIntT> dummy(img.Frame());
	dummy.Fill(0);
	RealT sfirst = -Log(scl.Covariance().Det());
	RealT lfirst = -Log(lcl.Covariance().Det());
	for(Array2dIter2C<UIntT, TFVectorC<RealT,3> > it(dummy, img); it; it++)
	{
		VectorC v(it.Data2()[0], it.Data2()[1], it.Data2()[2]);
		RealT Jlips = lfirst - lcl.MahalanobisDistance(v);
		RealT Jskin = sfirst - scl.MahalanobisDistance(v);
		RealT j = Jlips - Jskin;
		if(j > 0)
			it.Data1() = 1;
	}	
	ConnectedComponentsC<UIntT> connected ; 
	// Apply the algorithm 
	Tuple2C<ImageC<UIntT>, UIntT> result = connected.Apply(dummy);
	SegmentationC sm(result);
	SArray1dC<UIntT> sarea = sm.Areas().Copy();
	sarea.Sort();
	//Look for the label corresponding to the second largest area (after skin)
	IntT lar = sarea[1];
	UIntT label = 0;
	for(IntT i = 0; i < (IntT)sm.Areas().Size(); i++)
	{
		if((sm.Areas())[i] == (UIntT)lar)
			label = i;
	}
	//Insert all label pixels to be lips
	DListC<Index2dC> lind;
	for(Array2dIterC<UIntT> it(result.Data1()); it; it++)
	{
		if((*it)== label)
			lind.Append(it.Index().Copy());
	}
	return lind;
}

DListC<Index2dC> JCC(DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > &clusts,const ImageC<VectorC> &img)
{
	//Show off the clustered Image
	UIntT lab = 0;
	ImageC<UIntT> m_img(img.Frame(),0);
	for(DLIterC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > it(clusts); it; it++)
	{
		for(DLIterC<Tuple2C<VectorC,Index2dC> > it2((*it).Data2()); it2; it2++)
		{
			m_img[(*it2).Data2()] = lab;
		}
		lab++;
	}
	SegmentationC dummyseg(m_img, lab);
	if(!Save("@X: Test Segmentation Image",dummyseg.RandomImage()))
		cerr<<"Failed to show segmentation image"<<endl;

	clusts.MergeSort(mscomp);
	ImageC<IntT> im(img.Frame());
	im.Fill(0);
	MeanCovarianceC scl = clusts.Nth(0).Data1();
	MeanCovarianceC lcl = clusts.Nth(1).Data1();
	ImageC<UIntT> dummy(img.Frame());
	dummy.Fill(0);
	RealT sfirst = -Log(scl.Covariance().Det());
	RealT lfirst = -Log(lcl.Covariance().Det());
	for(Array2dIter2C<UIntT, VectorC> it(dummy, img); it; it++)
	{
		VectorC v(it.Data2().Copy());
		RealT Jlips = lfirst - lcl.MahalanobisDistance(v);
		RealT Jskin = sfirst - scl.MahalanobisDistance(v);
		RealT j = Jlips - Jskin;
		if(j > 0)
			it.Data1() = 1;
	}	
	ConnectedComponentsC<UIntT> connected ; 
	// Apply the algorithm 
	Tuple2C<ImageC<UIntT>, UIntT> result = connected.Apply(dummy);
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
	DListC<Index2dC> lind;
	for(Array2dIterC<UIntT> it(result.Data1()); it; it++)
	{
		if((*it)== label)
			lind.Append(it.Index().Copy());
	}
	return lind;	
}

ImageC<UIntT> GetBinaryLipImage(const DListC<Index2dC> &lip_lst, const ImageRectangleC &imrec)
{
	cout<<"Inside Binary Lip Image"<<endl;
	ImageC<UIntT> out(imrec,0);
	for(DLIterC<Index2dC> it(lip_lst); it; it++)
	{	
		//cout<<"P"<<endl;
		out[*it] = 255;
	}
	if(!Save("@X:LB Image Inside CF",out)) cerr<<"Could not save Lip Binary Image"<<endl;
	return out; 	
}

//ComputeSegmentationQuality () - Given two binary images L and G
//Compute the quality Q = L&&G / L||G
//Please note that ground truth binary condition is 0 and 255
//Lip image binary condition is 0 and 1
RealT ComputeSegmentationQuality(const ImageC<UIntT> &gimg, const ImageC<UIntT> &limg)
{
	RealT numerator = 0.0;
	RealT denominator = 0.0;
	for(Array2dIter2C<UIntT, UIntT> it(gimg, limg); it; it++)
	{
		if((it.Data1() == 255)&&(it.Data2() == 255))
			numerator++;
		if((it.Data1() == 255)||(it.Data2() == 255))	
			denominator++;
	}
	//cout<<"N = "<<numerator<<"\t D = "<<denominator<<endl;
	return (numerator/denominator);
}

//Function to compute the J Function suggested by URS
RealT ComputeJFunction(const MeanCovarianceC &mc, const VectorC vec)
{
	RealT output = -Log(mc.Covariance().Det()) - mc.MahalanobisDistance(vec);
	return output;
}

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

ImageC<UIntT> LabelImage(const ImageC<RealT> &im)
{
	ImageC<UIntT> res_im(im.Frame(),0);
	for(Array2dIter2C<RealT, UIntT> it(im, res_im); it; it++)
	{
		if(it.Data1() < 0)
			it.Data2() = 0; //Label Skin Class As Zero
		else
			it.Data2() = 255; //Label Lip Class As One
	}
	return res_im;
}
