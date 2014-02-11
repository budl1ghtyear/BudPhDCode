//File: ColourVectorExtractor.cc
//Input: Include a file containing some colour information
//Output: An output text file which contains a listing of all the individual colour vectors in the image
//Author: Bud Goswami
//Date: 23.01.09

#include <fstream>

#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter2.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/RealHSVValue.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Image/RealYUVValue.hh"
#include "Ravl/IO.hh"
#include "Ravl/Matrix.hh"
#include "Ravl/Option.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/StdMath.hh"
#include "Ravl/Vector.hh"
#include "Ravl/Vector2d.hh"
using namespace RavlImageN;

Array1dC<RealRGBValueC> GetRGBVals(const FilenameC &imgname);
Array1dC<Vector2dC> GetNormRGVals(const FilenameC &imgname);
Array1dC<RealHSVValueC> GetHSVVals(const FilenameC &imgname);
Array1dC<RealYUVValueC> GetYUVVals(const FilenameC &imgname);
Array1dC<VectorC> GetCIELabVals(const FilenameC &imgname);
Array1dC<VectorC> GetPseudoHueValsArray(const FilenameC &imgname);
int main(int nargs,char **argv) {
  OptionC opt(nargs,argv);
  FilenameC qimg = opt.String("i","in.ppm","Input Query Image File");
  FilenameC res = opt.String("o","out.txt","Output Colour Vectors");
  UIntT pix_type = opt.Int("p",0,"Pixel Output Type - 0=RealRGBValueC, 1=NormalisedRG, 2=RealHSVValueC, 3=RealYUVValueC, 4=CIELab");
  opt.Check();
  
  switch(pix_type)
  {
  case 0:
  {
 	Array1dC<RealRGBValueC> rgb = GetRGBVals(qimg);
 	ofstream ofs(res);
 	ofs << rgb;
 	ofs.close();
 	break;
  }
  case 1:
  {
	Array1dC<Vector2dC> img = GetNormRGVals(qimg);
	ofstream ofs(res);
	ofs << img;
	ofs.close();
	break;
  }
  case 2:
  {
	Array1dC<RealHSVValueC> img = GetHSVVals(qimg);
	ofstream ofs(res);
	ofs << img;
	ofs.close();
	break;
  }
  case 3:
  {
	Array1dC<RealYUVValueC> img = GetYUVVals(qimg);
	ofstream ofs(res);
	ofs << img;
	ofs.close();
	break;
  }
  case 4:
  {
	Array1dC<VectorC> img = GetCIELabVals(qimg);
	ofstream ofs(res);
	ofs << img;
	ofs.close();
	break; 	  
  }
  case 5:
  {
	Array1dC<VectorC> img = GetPseudoHueValsArray(qimg);
	ofstream ofs(res);
	ofs << img;
	ofs.close();
	break; 	  
  }
  }
  return 0;
}

Array1dC<RealRGBValueC> GetRGBVals(const FilenameC &imgname)
{
	ImageC<RealRGBValueC> image;
	if(!Load(imgname,image)) 
	{
	    cerr << "Failed to load file '" << imgname << "' \n";
	}
	Array1dC<RealRGBValueC> res(image.Frame().Area());
	IndexC ind = 0;
	for(Array2dIterC<RealRGBValueC> it(image); it; it++)
	{
		res[ind++] = (*it).Copy();
	}
	return res;  	
}

Array1dC<Vector2dC> GetNormRGVals(const FilenameC &imgname)
{
	ImageC<RealRGBValueC> image;
	if(!Load(imgname,image)) 
	{
	    cerr << "Failed to load file '" << imgname << "' \n";
	}
	Array1dC<Vector2dC> res(image.Frame().Area());
	Vector2dC zer(0.0,0.0);
	res.Fill(zer);
	IndexC ind = 0;
	for(Array2dIterC<RealRGBValueC> it(image); it; it++)
	{
		RealRGBValueC pix = (*it);
		RealT intensity = pix.Y()*3.0;
		res[ind++].Row() = pix.Red() / intensity;
		res[ind].Col() = pix.Green() / intensity;
	}
	return res;  	
}

Array1dC<RealHSVValueC> GetHSVVals(const FilenameC &imgname)
{
	ImageC<RealHSVValueC> image;
	if(!Load(imgname,image)) 
	{
	    cerr << "Failed to load file '" << imgname << "' \n";
	}
	Array1dC<RealHSVValueC> res(image.Frame().Area());
	IndexC ind = 0;
	for(Array2dIterC<RealHSVValueC> it(image); it; it++)
	{
		res[ind++] = (*it);
	}
	return res;  	
}

Array1dC<RealYUVValueC> GetYUVVals(const FilenameC &imgname)
{
	ImageC<RealYUVValueC> image;
	if(!Load(imgname,image)) 
	{
	    cerr << "Failed to load file '" << imgname << "' \n";
	}
	Array1dC<RealYUVValueC> res(image.Frame().Area());
	IndexC ind = 0;
	for(Array2dIterC<RealYUVValueC> it(image); it; it++)
	{
		res[ind++] = (*it);
	}
	return res;  	
}

Array1dC<VectorC> GetCIELabVals(const FilenameC &imgname)
{
	ImageC<RealRGBValueC> image;
	if(!Load(imgname,image)) 
	{
	    cerr << "Failed to load file '" << imgname << "' \n";
	}
	VectorC zer(0.0,0.0,0.0);
	Array1dC<VectorC> resrgb(image.Frame().Area());
	resrgb.Fill(zer);
	IndexC ind = 0;
	for(Array2dIterC<RealRGBValueC> imt(image); imt; imt++)
	{
		RealRGBValueC pix = (*imt);
		VectorC vec(pix.Red(), pix.Green(), pix.Blue());
		resrgb[ind++] = vec.Copy();
	}
	//Declare Standard Values
	VectorC refwhitergb(255.0,255.0,255.0);
	MatrixC xyzconv(0.490,0.310,0.200,0.177,0.813,0.011,0.000,0.010,0.990);
	//Obtain Reference White in XYZ
	VectorC refwhitexyz = xyzconv*refwhitergb;
	Array1dC<VectorC> rescielab(resrgb.Size());
	rescielab.Fill(zer);
	for(Array1dIter2C<VectorC, VectorC> it(resrgb, rescielab); it; it++)
	{
		//For each pixel, obtain the XYZ values
		VectorC xyzpix = xyzconv * it.Data1();
		VectorC cielabres(0.0,0.0,0.0);
		// Implement a loop to simultaneously calculate the X',Y',Z' values and obtain K1,K2,K3.
		for(IndexC i = 0; i < xyzpix.Size(); i++)
		{
			//Now xyzpix contains the normalised xyz values (X', Y', Z')
			xyzpix[i] /= refwhitexyz[i];
			if(i==1)
			{
				//Perform L* Calculations given Y'
				if(xyzpix[i] > 0.008856)
					cielabres[0] = (116.0 * Cbrt(xyzpix[i])) - 16.0;
				else
					cielabres[0] = 903.30 * xyzpix[i];
			}
			//Now Perform Kn calculations
			if(xyzpix[i] <= 0.008856)
				xyzpix[i] = (7.787*xyzpix[i]) + (16.0/116.0);
		}
		cielabres[1] = 500.0*(Cbrt(xyzpix[0]) - Cbrt(xyzpix[1])); //Calculate a*
		cielabres[2] = 200.0*(Cbrt(xyzpix[1]) - Cbrt(xyzpix[2])); //Calculate b*
		it.Data2() = cielabres.Copy();
	}
	return rescielab;
}

Array1dC<VectorC> GetPseudoHueValsArray(const FilenameC &imgname)
{
	//Load the image as an RGBValueC Image
	ImageC<RealRGBValueC> rgbimg;
	if(!Load(imgname, rgbimg)) cerr<<"Could not load RGB Image"<<endl;
	//Now use the PseudoHue Values method and populate an array with the observations
	Array1dC<VectorC> res(rgbimg.Frame().Area());
	UIntT i = 0;
	for(Array2dIterC<RealRGBValueC> it(rgbimg); it; it++)
	{
		RealT hue = (RealT)((*it).Red() / (RealT)((*it).Red()+(*it).Green()));
		VectorC vec(hue, (*it).Y());
		res[i++] = vec.Copy();
	}
	return res;
}
