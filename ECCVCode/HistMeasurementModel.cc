#include "HistMeasurementModel.hh"

RealT HistMeasurementModelC::Measure(ParticleC &pt)
{
	RealT weight = 0.;
	//BasisSplineC spl(ArrayToSArray_B(pt.GetState().GetObservationVector()));
	SArray1dC<Point2dC> curvepts = ArrayToSArray_B(pt.GetObservationVector());
	ImageRectangleC currentrect = GetImageRectangleFromPts(curvepts);
	//~ cout<<"Current Rect: "<<currentrect<<endl;
	Point3dC maskimg(8.,8.,8.);
	weight = BhattacharyaDistance(currentrect,imrect,maskimg);
	return weight;
}
RealT HistMeasurementModelC::MeasureQuality(SArray1dC<LinePP2dC> &norms)
{	
	RealT weight = 0.;
	SArray1dC<Point2dC> pts(norms.Size());
	for(SArray1dIter2C<LinePP2dC,Point2dC> it(norms,pts); it; it++)
		it.Data2() = it.Data1().FirstPoint().Copy();
	ImageRectangleC currentrect = GetImageRectangleFromPts(pts);
	Point3dC maskimg(8.,8.,8.);
	weight = BhattacharyaDistance(currentrect,imrect,maskimg);
	//~ cout<<"Weight = "<<weight<<endl;
	return weight;		
}
ImageRectangleC HistMeasurementModelC::GetImageRectangleFromPts(SArray1dC<Point2dC> const &pts)
{
	DListC<Point2dC> ppts;
	for(UIntT i = 0; i < pts.Size(); i++)
		ppts.Append(pts[i].Copy());
		
	Polygon2dC dummy(ppts);
	RealRange2dC rnge = dummy.BoundingRectangle();
	ImageRectangleC res(rnge.IndexRange());
	return res;
}

RealT HistMeasurementModelC::BhattacharyaDistance(ImageRectangleC const &rec1,ImageRectangleC const &rec2, Point3dC const &mask)
{
	//~ CvHistogram* hist1;
	//~ GetHistogram(im, rec1, hist1, mask);
	//~ CvHistogram* hist2;
	//~ GetHistogram(im, rec2, hist2, mask);
	//~ //CvHistogram *h1 = hist1;
	//~ //CvHistogram *h2 = hist2;
	//~ RealT res = cvCompareHist(hist1, hist2, CV_COMP_BHATTACHARYYA);
	//~ cvReleaseHist(&hist1);
	//~ cvReleaseHist(&hist2);
	//~ RealHistogram3dC hist1 = GetHistogram(im, rec1, mask);
	SArray3dC<UIntT> hist1 = GetHistogram(im, rec1, mask);
	//~ cout<<"Hist 1 Size: "<<hist1.Size()<<endl;
	//~ RealHistogram3dC hist2 = GetHistogram(im, rec2, mask);
	SArray3dC<UIntT> hist2 = GetHistogram(im, rec2, mask);
	//~ cout<<"Hist 1 Size: "<<hist2.Size()<<endl;
	RealT num(0.),den1(0.),den2(0.);
	RealT val(0.);
	for(UIntT i = 0; i < mask.X(); i++)
	{
		for(UIntT j = 0; j < mask.Y(); j++)
		{
			for(UIntT k = 0; k < mask.Z(); k++)
			{
				den1 += hist1[i][j][k];
				den2 += hist2[i][j][k];
				num = RealT((hist1[i][j][k])*(hist2[i][j][k]));
				//~ cout<<"Den 1 : "<<den1<<"\t Den2 : "<<den2<<"\t Num : "<<num<<endl;
				if(((den1*den2) == 0.)||(num == 0.))
					val+= 0.;
				else
					val += Sqrt(num) / Sqrt(den1 * den2);
			}
		}
	}
	RealT res = Sqrt(1 - val);	
	return res;	
}
