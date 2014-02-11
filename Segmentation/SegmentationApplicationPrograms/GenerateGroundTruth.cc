/*
#include "Ravl/Option.hh"
#include "Ravl/Image/ImagePointFeatureSet.hh"
#include "Ravl/HashIter.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/StrStream.hh"
#include "Ravl/XMLStream.hh"
#include "Ravl/EntryPnt.hh"
#include "Ravl/SArray1d.hh"
#include "Ravl/SArray1dIter.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/Index2d.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/DP/FileFormatIO.hh"
#include "Ravl/Polygon2d.hh"
using namespace RavlN;
using namespace RavlImageN;
*/ 

#include "XMLToGTConv.hh"
using namespace RavlN;
using namespace RavlImageN;

int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	StringC data = opt.String("i","test.xml","InputQueryXMLStream");
	//DirectoryC out_dir = opt.String("o","./Test/", "Location of stored Ground Truth");
	opt.Check();

	XMLToGTConv conv(data);

    /*
	XMLIStreamC strm(data);
	ImagePointFeatureSetC im_data;
	strm>>im_data;
	
	ImageC<RealRGBValueC> rgb_im;
	FilenameC f_name = "./"+im_data.ImageFile().NameComponent();
	if(!Load(f_name, rgb_im)) cerr<<"Could not load image file name"<<endl;
	if(!Save("@X:Input File", rgb_im)) cerr<<"Could not show image file name"<<endl;
	
	//Now lable lips in image
	RealRGBValueC lipcol(255,10,10);
	DListC<Point2dC> lip_pts;
	for( HashIterC<StringC,DListC<IntT> > it(im_data.SubsetIterator()); it; it++)
	{
		if(it.Key() == "Outer Lip 1")
		{
			for(DLIterC<IntT> dl(it.Data());dl; dl++ )
			{
				Point2dC pt = im_data.Position((*dl));
				lip_pts.Append(pt.Copy());
				Index2dC to_draw(pt.Row(), pt.Col());
				DrawCross(rgb_im,lipcol,to_draw,3);    
			}
		}
		else if(it.Key() == "Outer Lip 2")
		{
			for(DLIterC<IntT> dl(it.Data());dl; dl++ )
			{
				Point2dC pt = im_data.Position((*dl));
				lip_pts.Append(pt.Copy());
				Index2dC to_draw(pt.Row(), pt.Col());
				DrawCross(rgb_im,lipcol,to_draw,3);    
			}	
		}
	}
	if(!Save("@X:Lip File", rgb_im)) cerr<<"Could not show outer lip file name"<<endl;
	//Now do fill polygon stuff
	Polygon2dC lip_poly(lip_pts);
	ImageC<UIntT> binary(rgb_im.Frame(),0);
	for(Array2dIterC<UIntT> it(binary); it; it++)
	{
		if(lip_poly.Contains(it.Index()))
			(*it) = 255;
	}
	if(!Save("@X:Lip Binary File", binary)) cerr<<"Could not show binary file name"<<endl;
	*/
	return 0;
}
