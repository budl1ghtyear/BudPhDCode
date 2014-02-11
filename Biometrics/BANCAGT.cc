#include "BANCAGT.hh"

BANCAGT::BANCAGT(const FilenameC &fn)
{
	StringC out,in;
	IStreamC strm(fn);
	DListC<Point2dC> out_pts,in_pts;
	strm>>out>>out_pts;
	//cout<<"Outer Lip Coords "<<outer<<endl;
	if(!(strm.IsEndOfStream()))
	{
		//cout<<"Inner Coords Exist"<<endl;
		strm>>in;
		//cout<<"In - "<< in <<endl;
		if(in.contains("Inner",0)!= 0)
			strm>>in_pts;	
		else
			inner.Append(Point2dC(0.0,0.0));	
	}
	else
	{
		inner.Append(Point2dC(0.0,0.0));
	}
	strm.Close();		
	outer = RearrangePoints(out_pts);
	inner = RearrangePoints(in_pts);
}

BANCAGT::BANCAGT(const FilenameC &fn, const ImageRectangleC &imrec):f_name(fn),lip_bin(imrec)
{
	StringC out,in;
	IStreamC strm(fn);
	DListC<Point2dC> out_pts,in_pts;
	strm>>out>>out_pts;
	//cout<<"Outer Lip Coords "<<outer<<endl;
	if(!(strm.IsEndOfStream()))
	{
		//cout<<"Inner Coords Exist"<<endl;
		strm>>in;
		//cout<<"In - "<< in <<endl;
		if(in.contains("Inner",0)!= 0)
			strm>>in_pts;	
		else
			inner.Append(Point2dC(0.0,0.0));	
	}
	else
	{
		inner.Append(Point2dC(0.0,0.0));
	}
	strm.Close();
	//cout<<"Outer Lip Boundary Before Rearranging "<<out_pts<<endl;
	//cout<<"Inner Lip Boundary Before Rearranging "<<in_pts<<endl;
	outer = RearrangePoints(out_pts);
	inner = RearrangePoints(in_pts);
	//Now call the apply method
	//cout<<"out is "<<out<<"\n in is "<<in<<endl;
	//Now replace the co-ordinates to reflect the correct coordinates
	lip_bin = Apply();
}

ImageC<UIntT> BANCAGT::Apply(void)
{
	
	//cout<<"Outer Lip Boundary "<<outer<<endl;
	Polygon2dC lip_outer(outer);
	ImageC<UIntT> binary(lip_bin.Frame(),0);
	//Draw example image
	/*for(DLIterC<Point2dC> it(outer); it; it++)
	{
		DrawCross(binary,(UIntT)128,Index2dC((*it).Row(),(*it).Col()),3);		
	}
	*/
	//if(!Save("@X:Binary GT Image To Start",binary)) cerr<< "Could not show binary image" << endl;
	for(Array2dIterC<UIntT> it(binary); it; it++)
	{
		Index2dC this_ind = it.Index();
		Point2dC this_pt(this_ind.Row(),this_ind.Col());
		if(lip_outer.Contains(this_pt))
			(*it) = 255;
	}
	//if(!Save("@X:Binary GT Image",binary)) cerr<< "Could not show binary image" << endl;
	if((inner.Size() > 1))
	{
		//cout<<"Inner Lip Boundary "<<inner<<endl;
		Polygon2dC lip_inner(inner);
		for(Array2dIterC<UIntT> it(binary); it; it++)
		{
			Index2dC this_ind = it.Index();
			Point2dC this_pt(this_ind.Row(),this_ind.Col());
			if(lip_inner.Contains(this_pt)) 
				(*it) = 0;
		}
	}
	//if(!Save("@X:Binary Inner GT Image",binary)) cerr<< "Could not show binary image" << endl;
	return binary;
}

bool mscomp(const Point2dC &dat1, const Point2dC &dat2)
{
	return (dat1.Col() <= dat2.Col() ? true:false);
}
bool mscomp2(const Point2dC &dat1, const Point2dC &dat2)
{
	return (dat1.Col() >= dat2.Col() ? true:false);
}

DListC<Point2dC> BANCAGT::RearrangePoints(DListC<Point2dC> &pts)
{
	UIntT outer_dim = pts.Size() / 2;
	for(UIntT i = 0 ; i < outer_dim; i++)
	{
		Point2dC first(pts.Nth(i+outer_dim).Row(),pts.Nth(i).Row());
		Point2dC second(pts.Nth(i+outer_dim).Col(),pts.Nth(i).Col());
		pts.Nth(i) = first.Copy();
		pts.Nth(i+outer_dim) = second.Copy();
	}	
	//We need more operations now
	//Sort out the two lists into upper and lower coordinates
	DListC<Point2dC> upper, lower;
	Polygon2dC lip(pts); Point2dC center = lip.Centroid();
	for(DLIterC<Point2dC> it(pts); it; it++)
	{
		if((*it).Row() < center.Row())
			upper.Append((*it).Copy());
		else
			lower.Append((*it).Copy());
	}
	//Now sort the upper and lower points accordingly
	upper.MergeSort(mscomp);
	lower.MergeSort(mscomp2);
	for(DLIterC<Point2dC> it(lower); it; it++)
	{
		upper.Append((*it).Copy());
	}
	return upper;
}
