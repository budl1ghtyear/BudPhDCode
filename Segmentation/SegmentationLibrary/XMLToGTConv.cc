#include "XMLToGTConv.hh"

XMLToGTConv::XMLToGTConv(const FilenameC &fn, const ImageRectangleC &imrec):f_name(fn),lip_bin(imrec)
{
	XMLIStreamC strm(fn);
	strm >> im_data;
	strm.Close();
	//Now call the apply method
	lip_bin = Apply();
}
ImageC<UIntT> XMLToGTConv::Apply(void)
{
	
	//ImageC<RealRGBValueC> rgb_im;
	//FilenameC f_name = "./" + im_data.ImageFile().NameComponent();
	//if(!Load(f_name, rgb_im)) cerr<<"Could not load image file name"<<endl;
	#ifdef DEBUG
	if(!Save("@X:Input File", rgb_im)) cerr<<"Could not show image file name"<<endl;
	#endif
	//Now label lips in image
	//RealRGBValueC lipcol(255,10,10);
	DListC<Point2dC> lip_pts;
	for( HashIterC<StringC,DListC<IntT> > it(im_data.SubsetIterator()); it; it++)
	{
		if(it.Key() == "Outer Lip 1")
		{
			for(DLIterC<IntT> dl(it.Data());dl; dl++ )
			{
				Point2dC pt = im_data.Position((*dl));
				lip_pts.Append(pt.Copy());
				//Index2dC to_draw(pt.Row(), pt.Col());
				//DrawCross(rgb_im,lipcol,to_draw,3);    
			}
		}
		else if(it.Key() == "Outer Lip 2")
		{
			for(DLIterC<IntT> dl(it.Data());dl; dl++ )
			{
				Point2dC pt = im_data.Position((*dl));
				lip_pts.Append(pt.Copy());
				//Index2dC to_draw(pt.Row(), pt.Col());
				//DrawCross(rgb_im,lipcol,to_draw,3);    
			}	
		}
	}
	#ifdef DEBUG
	if(!Save("@X:Lip File", rgb_im)) cerr<<"Could not show outer lip file name"<<endl;
	#endif
	//Now do fill polygon stuff
	Polygon2dC lip_poly(lip_pts);
	ImageC<UIntT> binary(lip_bin.Frame(),0);
	for(Array2dIterC<UIntT> it(binary); it; it++)
	{
		if(lip_poly.Contains(it.Index()) && (binary.Frame().Contains(it.Index())))
			(*it) = 255;
	}
	#ifdef DEBUG
	if(!Save("@X:Lip Binary File", binary)) cerr<<"Could not show binary file name"<<endl;
	#endif
	return binary;
}


