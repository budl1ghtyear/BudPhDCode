#include "TrackingFunctions.hh"

Tuple2C<ImageC<RealT>,Affine2dC> ConvertImages(const ImageC<RealT> & image, const Point2dC &nle,const Point2dC &nre)
{
    ImageC<RealT> normalisedImage;
/*
//!Gabor
    RealT rows = 160;
    RealT cols = 128;
//for leye=28c,55r
    RealT rowFrac = 0.34375;
    RealT colFrac = 0.21875;
*/
  RealT rows = 142;
  RealT cols =  120;
  RealT rowFrac = 0.35;
  RealT colFrac = 0.25;

//for FRGC 3d Data
//        RealT rowFrac = 0.3;
 //   RealT colFrac = 0.1;
    UIntT gaussOrder = 7;
        
    //: Workout a simple geometric normalisation
    //==========================================
    Point2dC le(rows * rowFrac, (1.0-colFrac) * cols);
    Point2dC re(rows * rowFrac, cols * colFrac);
    ImageRectangleC outRect(Round(rows), Round(cols));
    Vector2dC d1 = le - re;
    Vector2dC d2 = nle - nre;
    RealT rot = d2.Angle() - d1.Angle();     
    RealT scale = d2.Modulus() / d1.Modulus();  
    Matrix2dC rotm = Matrix2dC(Cos(rot) * scale,-Sin(rot) * scale, Sin(rot) * scale,Cos(rot) * scale);      
    Point2dC cent= ((nle + nre)/2);
    Point2dC dcent = rotm * (le + re)/2;
    Vector2dC off = (cent - dcent);
    Affine2dC tr(rotm,off);
    
    //: Smooth and geometric normalise
    //=================================
    GaussConvolve2dC<RealT> smooth(gaussOrder);
    WarpAffineC<RealT> warp(outRect,tr);
    ImageRectangleC irec = warp.InputRectangle();
    irec = irec.Expand(gaussOrder/2 + 2);
    irec.ClipBy(image.Frame());
    normalisedImage = smooth.Apply(ImageC<RealT>(image,irec));

//  warp.SetOutputRectangle(IndexRange2dC(144,120));
    normalisedImage = warp.Apply(normalisedImage);
    
    //: Lets display the image
    //========================
	Tuple2C<ImageC<RealT>,Affine2dC> out(normalisedImage, tr);
   return out;
}


Array1dC<Point2dC> ConvertVecToPointArray(const VectorC &v)
{
	UIntT size = v.Size()/2;
	Array1dC<Point2dC> out(size);
	for(UIntT i = 0; i < size; i++)
	{	
		Point2dC pt(v[i + size], v[i]);
		out[i] = pt.Copy();
	}
	return out;
}

ImageC<RealT> ConvertRGBToRealImage(const ImageC<RealRGBValueC> &im)
{
	ImageC<RealT> res(im.Frame(),0.0);
	for(Array2dIter2C<RealRGBValueC,RealT> it(im, res); it; it++)
	{
		it.Data2() = it.Data1().Y();
	}
	return res;
}

Array1dC<Index2dC> GetEpsBoundary(const Array1dC<Index2dC> &ar, ImageC<RealT> &mc, ImageC<RealT> &sk)
{
	//Compute the centre point
	Index2dC ctr(0,0);
	IndexC size = 0;
	for(Array1dIterC<Index2dC> it(ar); it; it++)
	{
		ctr = ctr + (*it);
		size++;
	}
	ctr.Row() /= (IndexC)ar.Size();
	ctr.Col() /= (IndexC)ar.Size();
	//ctr /= (IndexC)ar.Size();	
	Array1dC<Index2dC> output(ar.Size());
	for(Array1dIter2C<Index2dC, Index2dC> it(ar, output); it; it++)
	{
		LinePP2dC line(ctr,it.Data1());
		Array1dC<RealT> left(22);
		left.Fill(0.0);
		Array1dC<RealT> right(22);
		right.Fill(0.0);
		IndexC arind = 0;
		//Create individual values of the MSdistance functions
		for(RealT i = 0.0; i <= 2.0; )
		{
			Point2dC loc = line.Point(i);
			if(i < 0.5)
			{
				right[arind] = 0.0;
				left[arind] = mc[(Index2dC)loc];
			}
			else if(i > 1.5)
			{
				left[arind] = 0.0;
				right[arind] = sk[(Index2dC)loc];
			}
			else
			{
				left[arind] = mc[(Index2dC)loc];
				right[arind] = sk[(Index2dC)loc];
			}
			i = i + 0.1;
			arind ++;
		}
		//Create an array of the epsilon term
		Array1dC<RealT> epsilon(5,15);
		epsilon.Fill(0.0);
		IndexC min=0;
		RealT minsum = 0.0;
		for(IndexC i = 5; i <= 15; i++)
		{
			RealT sumleft = 0.0;
			RealT sumright = 0.0;
			//sumleftb and right
			for(IndexC j = 0; j < 5; j++)
			{
				sumleft += left[i - j];
				sumright += right[i+j];
			}
			if(i == 5)
				minsum = sumleft + sumright;
			else
			{
				if(minsum >= (sumleft + sumright))
				{
					minsum = sumleft + sumright;
					min = i;
				}
			}
		}
		//now min contains the location of the point that has the minimum epsilon value
		Index2dC lipbndry = (Index2dC)line.Point((RealT)(min/10.0));
		it.Data2() = lipbndry.Copy();
	}
	return output;
}

//Method to get lip contour boundary points
//Generate ellipse
//Using ellipse, get B-Spline
//Using B_Spline normal, generate line searches using the binary image
//Return new points
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

//TESTING FUNCTIONS
void ComputeTest(const SArray1dC<Point2dC> &tres, BSplineC &spl, const ImageC<RealRGBValueC> &rgbimg, const UIntT &res)
{
	RealT white = 255.0;
	//This is when we get observations directly from a tres file
	ImageC<RealT> gt_img(rgbimg.Frame(),0.0);
	for(UIntT i = 0; i < (tres.Size()-1); i++)
	{
		Point2dC cur(tres[i].Col(),tres[i].Row());
		Point2dC nxt(tres[i+1].Col(),tres[i+1].Row());
		DrawLine(gt_img,white,cur,nxt);
	} 
	DrawLine(gt_img,white,Point2dC(tres[tres.Size()-1].Col(),tres[tres.Size()-1].Row()),Point2dC(tres[0].Col(),tres[0].Row()));
	if(!Save("@X:GT_Image", gt_img)) cerr<<"Could not save GT_IMG"<<endl;
	//Just for visual purposes, also draw on the same thing, the result found using crosses
	Array1dC<Point2dC> spl_pts = spl.RenderCurve(30);
	for(Array1dIterC<Point2dC> it(spl_pts); it; it++)
	{
		DrawCross(gt_img,white,Index2dC((*it).Row(), (*it).Col()),3);
	}
	if(!Save("@X:VIS_Image", gt_img)) cerr<<"Could not save GT_IMG"<<endl;
	//Now we have the two drawn curves, do a curve normal search for the ground-truth lip boundary
	/*RealT incr = 1.0/(RealT)res;	//curve resolution increment parameter
	//IntT search = 0;
	for(RealT pval = 0; pval < res; pval ++)		
	{
		RealT par = (RealT)pval*incr;
		LinePP2dC normal = bspl.GetNormalLine(par); //Get B-Spline Normal Line
		Point2dC start = normal.Point(0.0);
		//Now just search from the start to the end for the binary condition to end
		RealT index;
		//for(RealT i = -5; i < search; i++)
		IntT i = -5;
		while
		{
			Point2dC cur = normal.Point(i); Point2dC next = normal.Point(i+1);
			if(j_im.Frame().Contains(next))
			{
				if((j_im[cur] >= 0)&&(j_im[next] < 0))
				{
					index = i;
				}
			}
		}
		j_wt += GetGaussianValue(3,index);
	}*/
	
}

