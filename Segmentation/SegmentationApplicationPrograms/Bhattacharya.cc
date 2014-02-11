//File: Bhattacharya.cc
//Input: Include a directory containing the Lip Colour Example Images, and the type of pixel rep
//Output: A screen output of the diversity measure between each
//Author: Bud Goswami
//Date: 25.02.09

#include <fstream>

//DATA STRUCTURES
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter2.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Matrix.hh"
#include "Ravl/Vector.hh"
#include "Ravl/Vector2d.hh"
//IMAGE STUFF
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/RealHSVValue.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Image/RealYUVValue.hh"
#include "Ravl/IO.hh"
#include "Ravl/Option.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/StdMath.hh"
//PATTERN RECOGNITION STUFF
#include "Ravl/SumsNd2.hh"
#include "Ravl/MeanCovariance.hh"
#include "Ravl/MatrixRS.hh"
#include "Ravl/Math.hh"
#include "Ravl/PatternRec/DesignKMeans.hh"
#include "Ravl/PatternRec/DesignGaussianMixture.hh"
#include "Ravl/PatternRec/Sample.hh"
#include "Ravl/PatternRec/DistanceEuclidean.hh"
#include "Ravl/PatternRec/GaussianMixture.hh"
using namespace RavlImageN;
using namespace RavlConstN;

ImageC<RealRGBValueC> GetRGBVals(const FilenameC &imgname);
ImageC<VectorC> GetNormRGVals(const FilenameC &imgname);
ImageC<RealHSVValueC> GetHSVVals(const FilenameC &imgname);
ImageC<RealYUVValueC> GetYUVVals(const FilenameC &imgname);
ImageC<VectorC> GetCIELabVals(const FilenameC &imgname);
ImageC<VectorC> GetPseudoHueVals(const FilenameC &imgname);
RealT ComputeB(const MeanCovarianceC &skin, const MeanCovarianceC &lip);

template <typename T>
MeanCovarianceC ComputeMC(const ImageC<T> &img);
template <typename T>
RealT ComputePF(ImageC<T> &sk, ImageC<T> &lp, const UIntT &mixes);

RealT ComputeIntGauss(const IntT &dim, const MeanCovarianceC &mc_one, const MeanCovarianceC &mc_two, const RealT &m_one, const RealT &m_two);

RealT PolynomialSolver(const IntT &dim, const SArray1dC<MeanCovarianceC> &sk_mc, const SArray1dC<RealT> &sk_wt, const SArray1dC<MeanCovarianceC> &lp_mc, const SArray1dC<RealT> &lp_wt, const UIntT &mixes);

int main(int nargs,char **argv) 
{
  	OptionC opt(nargs,argv);
  	DirectoryC skin_dir = opt.String("s","Skin/","Skin Subjects");
	DirectoryC lip_dir = opt.String("l","Lip/","Lip Subjects");
	UIntT pix_type = opt.Int("p",0,"Pixel Output Type - 0=RealRGBValueC, 1=NormalisedRG, 2=RealHSVValueC, 3=RealYUVValueC, 4=CIELab, 5 = PseudoHue");
	UIntT nval = opt.Int("n",3,"Number of clusters you want to instantiate for GMM");
	opt.Check();
  
  
   //CREATE A LIST OF PPM FILES FOR EACH DIRECTORY
	DListC<StringC> skin_files = skin_dir.FiltList("*.ppm");
	DListC<StringC> lip_files = lip_dir.FiltList("*.ppm");
   
  	RealT bdist = 0.0;
  	for(DLIterC<StringC> d(skin_files); d; d++)
  	{
  		FilenameC sk_file = skin_dir + (*d).Copy();
  		FilenameC lp_file = lip_files.First().Copy();
  		for(DLIterC<StringC> it(lip_files); it; it++)
  		{
  			if((sk_file.NameComponent()[2] == (*it)[2])&&(sk_file.NameComponent()[3] == (*it)[3])&&(sk_file.NameComponent()[4] == (*it)[4]))
  			{
  				lp_file = lip_dir + (*it).Copy();
  				break;
  			}
  		}
  		//cout<<"Comparing "<<sk_file<<" with "<<lp_file<<endl; //Now we are comparing two same images
  		cout<<"File = "<<sk_file.BaseNameComponent()<<endl;
  		switch(pix_type)
  		{
  			case 0:
  			{
  				ImageC<RealRGBValueC> sk_img = GetRGBVals(sk_file);
  				ImageC<RealRGBValueC> lp_img = GetRGBVals(lp_file);
  				if(!Save("@X: Skin Image - ",sk_img)) cerr<<"Could not show skin image"<<endl;
  				if(!Save("@X: Lip Image - ",lp_img)) cerr<<"Could not show lip image"<<endl;
  				MeanCovarianceC skin = ComputeMC(sk_img);
  				MeanCovarianceC lip = ComputeMC(lp_img);
  				bdist = ComputeB(skin, lip);
  				RealT d = ComputePF(sk_img,lp_img,nval);
  				cout<<"Patrick-Fischer Distance = "<<d<<endl;
  				break;
  			}
 			case 1:
  			{
  				ImageC<VectorC> sk_img = GetNormRGVals(sk_file);
  				ImageC<VectorC> lp_img = GetNormRGVals(lp_file);
  				//if(!Save("@X: Skin Image - ",sk_img)) cerr<<"Could not show skin image"<<endl;
  				//if(!Save("@X: Lip Image - ",lp_img)) cerr<<"Could not show lip image"<<endl;
   			MeanCovarianceC skin = ComputeMC(sk_img);
  				MeanCovarianceC lip = ComputeMC(lp_img);
  				bdist = ComputeB(skin, lip); 	
  				RealT d = ComputePF(sk_img,lp_img,nval);		
  				cout<<"Patrick-Fischer Distance = "<<d<<endl;	
  				break;
  			}
  			case 2:
  			{
  				ImageC<RealHSVValueC> sk_img = GetHSVVals(sk_file);
  				ImageC<RealHSVValueC> lp_img = GetHSVVals(lp_file);
  				if(!Save("@X: Skin Image - ",sk_img)) cerr<<"Could not show skin image"<<endl;
  				if(!Save("@X: Lip Image - ",lp_img)) cerr<<"Could not show lip image"<<endl;
   			MeanCovarianceC skin = ComputeMC(sk_img);
  				MeanCovarianceC lip = ComputeMC(lp_img);
  				bdist = ComputeB(skin, lip); 				
  				RealT d = ComputePF(sk_img,lp_img,nval);
  				cout<<"Patrick-Fischer Distance = "<<d<<endl;
   			break;
  			}
  			case 3:
  			{
  				ImageC<RealYUVValueC> sk_img = GetYUVVals(sk_file);
  				ImageC<RealYUVValueC> lp_img = GetYUVVals(lp_file);
  				if(!Save("@X: Skin Image - ",sk_img)) cerr<<"Could not show skin image"<<endl;
  				if(!Save("@X: Lip Image - ",lp_img)) cerr<<"Could not show lip image"<<endl;
   			MeanCovarianceC skin = ComputeMC(sk_img);
  				MeanCovarianceC lip = ComputeMC(lp_img);
  				bdist = ComputeB(skin, lip); 				
  				RealT d = ComputePF(sk_img,lp_img,nval);
  				cout<<"Patrick-Fischer Distance = "<<d<<endl;
   			break;
  			}
 			case 4:
  			{
  				ImageC<VectorC> sk_img = GetCIELabVals(sk_file);
  				ImageC<VectorC> lp_img = GetCIELabVals(lp_file);
  				//if(!Save("@X: Skin Image - ",sk_img)) cerr<<"Could not show skin image"<<endl;
  				//if(!Save("@X: Lip Image - ",lp_img)) cerr<<"Could not show lip image"<<endl;
   			MeanCovarianceC skin = ComputeMC(sk_img);
  				MeanCovarianceC lip = ComputeMC(lp_img);
  				bdist = ComputeB(skin, lip); 				
  				RealT d = ComputePF(sk_img,lp_img,nval);
  				cout<<"Patrick-Fischer Distance = "<<d<<endl;
  				break;
  			}
 			case 5:
  			{
  				ImageC<VectorC> sk_img = GetPseudoHueVals(sk_file);
  				ImageC<VectorC> lp_img = GetPseudoHueVals(lp_file);
  				//if(!Save("@X: Skin Image - ",sk_img)) cerr<<"Could not show skin image"<<endl;
  				//if(!Save("@X: Lip Image - ",lp_img)) cerr<<"Could not show lip image"<<endl;
   			MeanCovarianceC skin = ComputeMC(sk_img);
  				MeanCovarianceC lip = ComputeMC(lp_img);
  				bdist = ComputeB(skin, lip); 				
  				RealT d = ComputePF(sk_img,lp_img,nval);
  				cout<<"Patrick-Fischer Distance = "<<d<<endl;
  				break;
  			}
 			
    		
  		}
  		//cout<<"Computed Bhattacharya Distance = "<<bdist<<endl;
  	}
   return 0;
}

ImageC<RealRGBValueC> GetRGBVals(const FilenameC &imgname)
{
	ImageC<RealRGBValueC> image;
	if(!Load(imgname,image)) 
	{
	    cerr << "Failed to load file '" << imgname << "' \n";
	}
	return image;  	
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

ImageC<RealHSVValueC> GetHSVVals(const FilenameC &imgname)
{
	ImageC<RealHSVValueC> image;
	if(!Load(imgname,image)) 
	{
	    cerr << "Failed to load file '" << imgname << "' \n";
	}
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

ImageC<VectorC> GetCIELabVals(const FilenameC &imgname)
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
template <typename T>
MeanCovarianceC ComputeMC(const ImageC<T> &img)
{
	//Given an image, treat the pixels as having a Vector interface
	//Compute the MeanCovariance and return
	bool s_stats = true; //for SumsN2dC object, normalise sample stats for comparison
	UIntT dim = img[img.Frame().TopLeft()].Size();
	SumsNd2C sum(dim);
	for(Array2dIterC<T> it(img); it; it++)
	{
		VectorC vec(dim);
		for(UIntT i = 0; i < dim; i++)
		{
			vec[i] = (*it)[i];
		}
		sum += vec.Copy();
	}
	return MeanCovarianceC(sum.MeanCovariance(s_stats)).Copy();
}

RealT ComputeB(const MeanCovarianceC &skin, const MeanCovarianceC &lip)
{
	MatrixRSC new_cov = MatrixRSC(skin.Covariance() + lip.Covariance())*(RealT)(0.5);
	MeanCovarianceC new_mc((skin.Number()+lip.Number()),lip.Mean(), new_cov);
	//cout<<"NEW mC "<<new_mc<<endl;
	RealT out = Log(new_cov.Det()/Sqr(skin.Covariance().Det() + lip.Covariance().Det()));
	//cout<<"ln term = "<<out<<endl;
	out += new_mc.MahalanobisDistance(skin.Mean());
	return out;
}
template <typename T>
RealT ComputePF(ImageC<T> &sk,ImageC<T> &lp, const UIntT &mixes)
{
	//compute patrick fischer distance
	//create 3 data structures to store skin, lip and combine data
	SampleC<VectorC> skin(sk.Size()), lip(lp.Size());
	for(Array2dIterC<T> t(sk); t; t++)
	{
		VectorC v((*t).Size());
		for(UIntT j = 0; j < v.Size(); j++)
		{
			v[j] = (*t)[j];
		}
		skin.Append(v.Copy());
	}
	bool isDiag = false;
	DistanceEuclideanC dist;
	DesignKMeansC km(mixes,dist);
	DesignGaussianMixtureC sk_gm(mixes, isDiag,km);
	GaussianMixtureC sk_fn = sk_gm.Apply(skin);
	for(Array2dIterC<T> it(lp); it; it++)
	{
		VectorC v((*it).Size());
		for(UIntT j = 0; j < v.Size(); j++)
		{
			v[j] = (*it)[j];
		}
		lip.Append(v.Copy());
	}
	DesignKMeansC lpkm(mixes,dist);
	DesignGaussianMixtureC lp_gm(mixes, isDiag,km);
	GaussianMixtureC lp_fn = lp_gm.Apply(lip);
	
	RealT ans = 0.0;
	UIntT dim = sk[sk.Frame().TopLeft()].Size();
	SArray1dC<MeanCovarianceC> sk_mc = sk_fn.MeanCovariances();
	SArray1dC<RealT> sk_wt = sk_fn.Weights();
	SArray1dC<MeanCovarianceC> lp_mc = lp_fn.MeanCovariances();
	SArray1dC<RealT> lp_wt = lp_fn.Weights();
	ans = PolynomialSolver(dim, sk_mc, sk_wt, lp_mc, lp_wt, mixes);
	return Sqrt(ans);
}

RealT ComputeIntGauss(const IntT &dim, const MeanCovarianceC &mc_one, const MeanCovarianceC &mc_two, const RealT &m_one, const RealT &m_two)
{
	RealT res = 0.0;
	//Case if the two Gaussians are still the same:
	if((mc_one.Mean() == mc_two.Mean())&&(mc_one.Covariance() == mc_two.Covariance()))
	{
		MatrixRSC cov(mc_one.Covariance() + mc_two.Covariance());
		//Res = {{(2Pi)^d}*|2Sigma|}^-1/2
		res = (m_one*m_two)* RealT(Pow((Pow((2.0*pi),dim)* cov.Det()),-0.50));
	}
	else
	{
		MatrixRSC cov(mc_one.Covariance() + mc_two.Covariance());
		//Res = {{(2Pi)^d}*|2Sigma|}^-1/2 * exp{-0.50(m2 - m1)^T*Sigma^-1*(m2 - m1)}
		res = (m_one * m_two)*(Pow((Pow((2.0*pi),dim)* cov.Det()),-0.50)); 
		MeanCovarianceC dmc((mc_one.Number() + mc_two.Number()),mc_one.Mean(),cov);
		RealT sec = Exp(-0.50 * dmc.MahalanobisDistance(mc_two.Mean()));
		res *= sec;			
	}
	return res;
}

//We have an eqn of the form ((x1 + x2 +...+xN)-(y0 + y1 + ... + yN))^2
//Solution - {x1(f(x))+x2(f(x))+...+xN(f(x))} - 2{x1(f(y))+x2(f(y))+...+xN(f(y)} + {y0(f(y)) + y1(f(y)) + ... +yN(f(y))}
RealT PolynomialSolver(const IntT &dim, const SArray1dC<MeanCovarianceC> &sk_mc, const SArray1dC<RealT> &sk_wt, const SArray1dC<MeanCovarianceC> &lp_mc, const SArray1dC<RealT> &lp_wt, const UIntT &mixes)
{
	RealT xterms = 0.0,midterms = 0.0,yterms = 0.0;
	//Terms
	for(UIntT i = 0; i < mixes; i++)
	{
		for(UIntT j = 0; j < mixes; j++)
		{
			xterms += ComputeIntGauss(dim, sk_mc[i],sk_mc[j],sk_wt[i],sk_wt[j]);
			midterms += ComputeIntGauss(dim, sk_mc[i],lp_mc[j],sk_wt[i],lp_wt[j]);
			yterms += ComputeIntGauss(dim, lp_mc[i],lp_mc[j],lp_wt[i],lp_wt[j]);
		}
	}
	cout<<"xterms = "<<xterms<<"\t 2.0 * midterms = "<<(2.0 * midterms)<<"\t y terms = "<<yterms<<endl;
	//return (xterms - (2.0 * midterms) + yterms);
	//New Patrick Fisher Distance Scaling Term
	return (RealT)((RealT)Sqrt((xterms - (2.0 * midterms) + yterms))/(RealT)Sqrt((xterms + yterms)));
}
