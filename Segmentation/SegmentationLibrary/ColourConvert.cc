//ColourConvert.cc
//Author - Bud Goswami
//Date - 18/03/09
//Contains the function definitions for the colour conversion functions

#include "ColourConvert.hh"

ImageC<RealRGBValueC> GetRGBVals(const FilenameC &imgname)
{
	ImageC<RealRGBValueC> image;
	if(!Load(imgname,image)) 
	{
	    cerr << "Failed to load file '" << imgname << "' \n";
	}
	return image;  	
}
ImageC<RealRGBValueC> GetRGBVals(const ImageC<RealRGBValueC> &img)
{
	return img.Copy();  	
}
ImageC<VectorC> GetNormRGVals(const FilenameC &imgname)
{
	ImageC<RealRGBValueC> image;
	if(!Load(imgname,image)) 
	{
	    cerr << "Failed to load file '" << imgname << "' \n";
	}
	ImageC<VectorC> res(image.Frame());
	VectorC zer(0.0,0.0);
	res.Fill(zer);
	for(Array2dIter2C<RealRGBValueC, VectorC> it(image, res); it; it++)
	{
		RealRGBValueC pix = it.Data1().Copy();
		RealT intensity = pix.Y()*3.0;
		VectorC vec( pix.Red() / intensity,pix.Green() / intensity);
		it.Data2() = vec.Copy();
	}
	return res;  	
}

ImageC<VectorC> GetNormRGVals(const ImageC<RealRGBValueC> &img)
{
	ImageC<VectorC> res(img.Frame());
	VectorC zer(0.0,0.0);
	res.Fill(zer);
	for(Array2dIter2C<RealRGBValueC, VectorC> it(img, res); it; it++)
	{
		RealRGBValueC pix = it.Data1().Copy();
		RealT intensity = pix.Y()*3.0;
		VectorC vec( pix.Red() / intensity,pix.Green() / intensity);
		it.Data2() = vec.Copy();
	}
	return res;  	
}

ImageC<RealHSVValueC> GetHSVVals(const FilenameC &imgname)
{
	ImageC<RealHSVValueC> image;
	if(!Load(imgname,image)) 
	{
	    cerr << "Failed to load file '" << imgname << "' \n";
	}
	return image;  	
}
ImageC<RealHSVValueC> GetHSVVals(const ImageC<RealRGBValueC> &im)
{
	ImageC<RealHSVValueC> image = RealRGBImageCT2RealHSVImageCT(im).Copy()  ;
	return image;  	
}

ImageC<RealYUVValueC> GetYUVVals(const FilenameC &imgname)
{
	ImageC<RealYUVValueC> image;
	if(!Load(imgname,image)) 
	{
	    cerr << "Failed to load file '" << imgname << "' \n";
	}
	return image;  	
}

ImageC<RealYUVValueC> GetYUVVals(const ImageC<RealRGBValueC> &im)
{
	ImageC<RealYUVValueC> image = RealRGBImageCT2RealYUVImageCT(im)    ;
	return image;  	
}

ImageC<VectorC> GetCIELabVals(const FilenameC &imgname)
{
	ImageC<RealRGBValueC> image;
	if(!Load(imgname,image)) 
	{
	    cerr << "Failed to load file '" << imgname << "' \n";
	}
	return GetCIELabVals(image);
}

ImageC<VectorC> GetCIELabVals(const ImageC<RealRGBValueC> &im)
{
	ImageC<RealRGBValueC> image = im.Copy();
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
	ImageC<VectorC> out(image.Frame());
	ind = 0;
	for(Array2dIterC<VectorC> it(out); it; it++)
	{
		(*it) = rescielab[ind++].Copy();
	}
	return out;
}


ImageC<VectorC> GetPseudoHueVals(const FilenameC &imgname)
{
	ImageC<RealRGBValueC> rgbimage;
	if(!Load(imgname, rgbimage)) cerr<<"Could not load "<<imgname<<endl;
	return GetPseudoHueVals(rgbimage);
}
ImageC<VectorC> GetPseudoHueVals(const ImageC<RealRGBValueC> &im)
{
	ImageC<RealRGBValueC> rgbimage = im.Copy();
	//Given RGB Image, perform pseudohue and luminance calculations and store in 2D vector
	ImageC<VectorC> out(rgbimage.Frame());
	for(Array2dIter2C<RealRGBValueC, VectorC> it(rgbimage,out); it; it++)
	{
		RealT hue = (RealT)(it.Data1().Red() / (RealT)(it.Data1().Red()+it.Data1().Green()));
		VectorC v(hue, it.Data1().Y());
		it.Data2() = v.Copy();
	}
	return out;
}
