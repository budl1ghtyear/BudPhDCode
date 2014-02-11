#include "Align.hh"

DListC<SArray1dC<Point2dC> > AlignC::PerformSimilarityAlignment(const DListC<SArray1dC<Point2dC> > &pts, const SArray1dC<Point2dC> &ref)
{
	//Perform Similarity Projection of all the points
	DListC<SArray1dC<Point2dC> > res;
	for(DLIterC<SArray1dC<Point2dC> > it(pts); it; it++)
	{
		//Compute the similarity projection and project the points
		bool forceUnitScale = false;
		Affine2dC transform;
		FitSimilarity(ref,(*it),transform,forceUnitScale); //Find the affine transform that projects point2 to point1   		
		SArray1dC<Point2dC> res_pts = AffineProject((*it), transform).Copy();
		res.Append(res_pts.Copy());
	}
	return res;
}
DListC<SArray1dC<Point2dC> > AlignC::PerformAffineAlignment(const DListC<SArray1dC<Point2dC> > &pts, const SArray1dC<Point2dC> &ref)
{
	DListC<SArray1dC<Point2dC> > res;
	for(DLIterC<SArray1dC<Point2dC> > it(pts); it; it++)
	{
		//Compute the similarity projection and project the points
		Affine2dC transform = 
		FitAffine((*it),ref); //Find the affine transform that projects point2 to point1   		
		SArray1dC<Point2dC> res_pts = AffineProject((*it), transform).Copy();
		res.Append(res_pts.Copy());
	}
	return res;	
}
DListC<SArray1dC<Point2dC> > AlignC::ComputeDeviation(const DListC<SArray1dC<Point2dC> > &pts,  const SArray1dC<Point2dC> &ref)
{
	DListC<SArray1dC<Point2dC> > res;
	for(DLIterC<SArray1dC<Point2dC> > it(pts); it; it++)
	{
		res.Append((SArray1dC<Point2dC>)((*it) - ref).Copy());
	}
	return res;	
}

SArray1dC<Point2dC> AlignC::AffineProject(const SArray1dC<Point2dC> &pts, const Affine2dC &aff)
{
	SArray1dC<Point2dC> res(pts.Size());
	for(UIntT i = 0; i < pts.Size(); i++)
	{
		Vector2dC vec(pts[i].Row(),pts[i].Col());
		Vector2dC proj_vec = aff * vec;
		Point2dC proj_pt(proj_vec.Row(),proj_vec.Col());
		res[i] = proj_pt.Copy();
	}
	return res;	
}
