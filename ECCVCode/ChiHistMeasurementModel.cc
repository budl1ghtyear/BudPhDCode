#include "ChiHistMeasurementModel.hh"

RealT ChiHistMeasurementModelC::Measure(ParticleC &pt)
{
	RealT weight = 0.;
	BasisSplineC spl(ArrayToSArray_B(pt.GetState().GetObservationVector()));
	SArray1dC<Point2dC> curvepts = spl.GetCurvePoints();
	ImageRectangleC currentrect = GetImageRectangleFromPts(curvepts);
	Point3dC maskimg(8.,8.,8.);
	weight = ChiDistance(currentrect,imrect,maskimg);
	return weight;
}
RealT ChiHistMeasurementModelC::MeasureQuality(SArray1dC<LinePP2dC> &norms)
{	
	RealT weight = 0.;
	SArray1dC<Point2dC> pts(norms.Size());
	for(SArray1dIter2C<LinePP2dC,Point2dC> it(norms,pts); it; it++)
		it.Data2() = it.Data1().FirstPoint().Copy();
	ImageRectangleC currentrect = GetImageRectangleFromPts(pts);
	Point3dC maskimg(8.,8.,8.);
	weight = ChiDistance(currentrect,imrect,maskimg);
	return weight;		
}
ImageRectangleC ChiHistMeasurementModelC::GetImageRectangleFromPts(SArray1dC<Point2dC> const &pts)
{
	DListC<Point2dC> ppts;
	for(UIntT i = 0; i < pts.Size(); i++)
		ppts.Append(pts[i].Copy());
		
	Polygon2dC dummy(ppts);
	RealRange2dC rnge = dummy.BoundingRectangle();
	ImageRectangleC res(rnge.IndexRange());
	return res;
}

RealT ChiHistMeasurementModelC::ChiDistance(ImageRectangleC const &rec1,ImageRectangleC const &rec2, Point3dC const &mask)
{
	SArray3dC<UIntT> hist1 = GetHistogram(im, rec1, mask);
	SArray3dC<UIntT> hist2 = GetHistogram(im, rec2, mask);
	RealT num(0.),den(0.);
	RealT val(0.);
	for(UIntT i = 0; i < mask.X(); i++)
	{
		for(UIntT j = 0; j < mask.Y(); j++)
		{
			for(UIntT k = 0; k < mask.Z(); k++)
			{
				den = hist1[i][j][k] + hist2[i][j][k];
				num = Sqr(RealT((hist1[i][j][k])-(hist2[i][j][k])));
				if((den) == 0.)
					val+= 0.;
				else
					val += num / den;
			}
		}
	}
	RealT res = Sqrt(1 - val);	
	return res;	
}
