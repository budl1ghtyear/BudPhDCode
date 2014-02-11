#include "ProcrustesAlign.hh"
ProcrustesAlignC::ProcrustesAlignC(const Array1dC<Point2dC> &ref)
{
	//perform data structure change:
	SArray1dC<Point2dC> dummy(ref.Size());
	for(UIntT i = 0; i < ref.Size(); i++)
	{
		dummy[i] = ref[i].Copy();
	}
	ref_shape = dummy.Copy();
}
ProcrustesAlignC::ProcrustesAlignC(const SArray1dC<Point2dC> &ref)
{
	ref_shape = ref.Copy();
}
ProcrustesAlignC::ProcrustesAlignC(const VectorC &ref)
{
	//Remember the data format inside a vector:{X0....Xm,Y0...Ym} i.e. {Col[0]....Col[m],Row[0]...Row[m]}
	UIntT dim = (UIntT)(ref.Size() / 2);
	SArray1dC<Point2dC> dummy(dim);
	for(UIntT i = 0; i < dummy.Size(); i++)
	{
		Point2dC pt(ref[i+dim],ref[i]);
		dummy[i] = pt.Copy();
	}
	ref_shape = dummy.Copy();	
}
//Apply method:
DListC<SArray1dC<Point2dC> > ProcrustesAlignC::Apply(const DListC<SArray1dC<Point2dC> > &pts)
{
	//Assign a new mean shape
	UIntT iter_count = 0;
	SArray1dC<Point2dC> newref_shape = ComputeMeanShape(pts).Copy();
	DListC<SArray1dC<Point2dC> > aligned_pts = pts.Copy();
	do
	{
		for(DLIterC<SArray1dC<Point2dC> > it(aligned_pts); it; it++)
		{
			SArray1dC<Point2dC> aligned_point = AlignShape(*it).Copy(); //Compute the aligned point
			(*it) = aligned_point; // assign the point to the aligned point
		}
		ref_shape = newref_shape;
		newref_shape = ComputeMeanShape(aligned_pts).Copy();
		cout<<"Iteration Count = "<<iter_count++<<"\n Ref_shape = "<<ref_shape<<endl;
	}while(ref_shape != newref_shape);
	return aligned_pts;
}

SArray1dC<Point2dC> ProcrustesAlignC::AlignShape(const SArray1dC<Point2dC> &pts)
{
	bool forceUnitScale = false;
	Affine2dC transform;
	FitSimilarity(ref_shape,pts,transform,forceUnitScale); //Find the affine transform that projects point2 to point1   	
	SArray1dC<Point2dC> res = AffineProject(pts, transform).Copy();
	return res;
}

SArray1dC<Point2dC> ProcrustesAlignC::AffineProject(const SArray1dC<Point2dC> &pts, const Affine2dC &aff)
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

SArray1dC<Point2dC> ProcrustesAlignC::ComputeMeanShape(const DListC<SArray1dC<Point2dC> > &pts)
{
	SArray1dC<Point2dC> tot(pts.First().Size());
	tot.Fill(Point2dC(0.0,0.0));
	for(DLIterC<SArray1dC<Point2dC> > it(pts);it;it++)
	{
		tot += (*it);
	}
	Point2dC elem((RealT)(pts.Size()),(RealT)(pts.Size()));
	tot /= elem;
	return (tot);
}		
