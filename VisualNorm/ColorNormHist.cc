///////////////////////////////////////////////////////////////////
//! file="ColorNormHist.cc"
//! author="Xuan Zou"

#include         "ColorNorm.hh"
#include 	"Ravl/Vector3d.hh"
#include 	"Ravl/Array2dIter.hh"
#include 	"Ravl/Array2dIter2.hh"
#include 	"Ravl/Array1dIter2.hh"
#include 	"Ravl/LeastSquares.hh"
#include 	"Ravl/Array1dIter.hh"
#include 	"Ravl/SArray1dIter.hh"

using namespace RavlN;
using namespace RavlImageN;

bool ColorNormHistStatC::GetHistStatistics(const ImageC<ByteRGBValueC> &im, Vector3dC &meanVec,  Vector3dC &stdVec)
{	RealT meanR =0, meanG =0, meanB=0;
	RealT stdR=0, stdG=0, stdB=0;
	for(Array2dIterC<ByteRGBValueC> it(im); it; it++)
	{	meanR+=(*it).Red();	meanG+=(*it).Green();	meanB+=(*it).Blue();	
	}
	meanR/=im.Size(); meanG/=im.Size(); meanB/=im.Size();
	for(Array2dIterC<ByteRGBValueC> it(im); it; it++)
	{	stdR+=((*it).Red() - meanR)*((*it).Red() - meanR);
		stdG+=((*it).Green() - meanG)*((*it).Green() - meanG);
		stdB+=((*it).Blue() - meanB)*((*it).Blue() - meanB);
	}
	stdR = sqrt(stdR/im.Size());  stdG = sqrt(stdG/im.Size());  stdB = sqrt(stdB/im.Size());  
	meanVec=Vector3dC(meanR, meanG, meanB); stdVec=Vector3dC(stdR,stdG,stdB);
	return 1;
};

bool ColorNormHistStatC::GetHistStatistics(const ImageC<ByteRGBValueC> &im, const ImageRectangleC &imRec, Vector3dC &meanVec,  Vector3dC &stdVec)
{	ImageC<ByteRGBValueC> imCrop(im, imRec);
	ImageC<ByteRGBValueC> imRegion(imCrop.Frame());
	for(	Array2dIter2C<ByteRGBValueC, ByteRGBValueC> it(imCrop, imRegion); it;it++)
	{	it.Data2() = it.Data1(); }
	GetHistStatistics(imRegion, meanVec, stdVec);
	return 1;
}


ImageC<ByteRGBValueC> ColorNormHistStatC::Apply(const ImageC<ByteRGBValueC> &imRGB, const Vector3dC &meanVecDes, const Vector3dC & stdVecDes, const Vector3dC &meanVecSrc, const Vector3dC& stdVecSrc)
{	
	Vector3dC factorVec = stdVecDes/ stdVecSrc ;
	ImageC<ByteRGBValueC> imout (imRGB.Frame());
	
	for(Array2dIter2C<ByteRGBValueC, ByteRGBValueC> it(imout, imRGB) ; it;it++)
	{	RealT tempR = ( (RealT) it.Data2().Red() - meanVecSrc.X() )* factorVec.X() + meanVecDes.X() ;
		if (tempR>255) it.Data1().Red()=255; 
		else if(tempR<0) it.Data1().Red() =0;
		else it.Data1().Red()=(ByteT) tempR;
	
		RealT tempG =  ((RealT) it.Data2().Green() - meanVecSrc.Y() )* factorVec.Y() + meanVecDes.Y() ;
		if (tempG>255) it.Data1().Green()=255; 
		else if(tempG <0) it.Data1().Green() =0;
		else it.Data1().Green()=(ByteT) tempG;

		RealT tempB = ( (RealT) it.Data2().Blue() - meanVecSrc.Z())* factorVec.Z() + meanVecDes.Z() ;
		if (tempB>255) it.Data1().Blue() =255; 
		else if(tempB <0) it.Data1().Blue() =0;
		else it.Data1().Blue()=(ByteT) tempB;
	}
	return imout;
};

bool ColorNormHistShapeC::GetColorMap(const ImageC<ByteT> &imSrc, const ImageC<ByteT> &imDes, Array1dC<ByteT> &m_ColorMap)
{	Array1dC<ByteT> ColorMap(256);
	ColorMap.Fill(0);
	Array1dC<RealT>	histSrc(256),histDes(256);
	histSrc.Fill(0); histDes.Fill(0);
	for(Array2dIterC<ByteT> it(imSrc); it;it++)
		histSrc[*it] ++;
	for(Array2dIterC<ByteT> it(imDes); it;it++)
		histDes[*it] ++;
	histSrc =histSrc/imSrc.Size();
	histDes =histDes/imDes.Size();
	Array1dC<RealT> cum_Src(256), cum_Des(256);
	cum_Src.Fill(0); cum_Des.Fill(0);
	RealT temp=0;
	//cout<<"-------------------\n hist1"<<endl;
	for(Array1dIter2C<RealT, RealT> it(histSrc, cum_Src); it;it++)
	{	temp=temp+it.Data1(); it.Data2() =temp;
		//cout<<it.Data1()<<' ';
	}
	//cout<<"\n hist2\n";
	temp=0;
	for(Array1dIter2C<RealT, RealT> it(histDes, cum_Des); it;it++)
	{	temp=temp+it.Data1(); it.Data2() =temp;
		//cout<<it.Data1()<<' ';
	}
	//cout<<endl;
	IndexC IndMinNZ_s =0, IndMinNZ_d=0;
	RealT thresh_s= 1.0/imSrc.Size();
	//cout<<cum_Src[IndMinNZ_s] <<' '<<thresh_s<<' '<<cum_Src[IndMinNZ_s] -thresh_s<<endl;
		
	while (cum_Src[IndMinNZ_s] < thresh_s)
		IndMinNZ_s ++;
	RealT thresh_d= 1.0/imDes.Size();
	while (cum_Des[IndMinNZ_d] < thresh_d)
	    IndMinNZ_d ++;
	//cout<<"IndMinNZ_s="<<IndMinNZ_s<<" IndMinNZ_d="<<IndMinNZ_d<<endl;
	
	IndexC IndMax_s =IndMinNZ_s;
	while (cum_Src[IndMax_s]<1.0-thresh_s/2 && IndMax_s<254)
	{		//cout<<"IndMax_S = "<<IndMax_s<<" cum_s[IndMax_s] = "<<cum_Src[IndMax_s]<< " cum_s[IndMax_s]<1? "<<(cum_Src[IndMax_s]<1)<<" diff:"<<cum_Des[IndMax_s]-1.0<<endl;
			IndMax_s++;
	}
	//cout<<"IndMax_s: "<<IndMax_s<<endl;
	
	IndexC IndMax_d =IndMinNZ_d;    
	while (cum_Des[IndMax_d]<1.0-thresh_d/2 && IndMax_d<254)
	{	//cout<<"IndMax_d = "<<IndMax_d<<" cum_d[IndMax_d] = "<<cum_Des[IndMax_d]<< " cum_Des[IndMax_d]<1? "<<(cum_Des[IndMax_d]<1)<<" diff:"<<cum_Des[IndMax_d]-1.0<<endl;
		IndMax_d++;
	}
	//cout<<"IndMax_d "<<IndMax_d<<endl;
	
	ColorMap[IndMinNZ_s] =(ByteT) IndMinNZ_d;
	IndexC preInd_d =IndMinNZ_d; 
	for (IndexC ind_s = IndMinNZ_s+1; ind_s<=IndMax_s ;ind_s++)
	{	//cout<<"preInd_d="<<preInd_d<<endl;
		ColorMap[ind_s] =preInd_d;
		if ( preInd_d==IndMax_d)
		 {	//cout<<"ColorMap["<<ind_s<<"]="<<(UIntT)ColorMap[ind_s]<<endl;
		 	continue;	}
		IndexC     ind_d = preInd_d;
		while (ind_d < IndMax_d)
		{      if (cum_Src[ind_s] <= cum_Des[ind_d] || ( cum_Src[ind_s] >= cum_Des[ind_d] && cum_Src[ind_s] <= cum_Des[ind_d+1])  )
           		{	 ColorMap[ind_s] = ind_d; 
				preInd_d =ind_d; 
				//cout<<"ColorMap["<<ind_s<<"]="<<(UIntT)ColorMap[ind_s]<<endl;
				break;
			}
		        else
             			ind_d = ind_d +1;
		}
        }
	ColorMap[IndMax_s] = IndMax_d;
	//cout<<"mid portion done"<<endl;
	if(IndMinNZ_s>0)
	{	//RealT colorStep= (RealT)ColorMap[IndMinNZ_s]/(RealT)IndMinNZ_s;
		for (IndexC ind_s =0; ind_s<IndMinNZ_s ;ind_s++)
			ColorMap[ind_s] =ColorMap[IndMinNZ_s];//(ByteT)(colorStep*ind_s);
	}
	//cout<<"first portion interpolated"<<endl;
	//cout<<"ColorMap[IndMax_s]="<<(UIntT)ColorMap[IndMax_s]<<endl;
	if(IndMax_s <255)
	{	//RealT colorStep= (255-(RealT)ColorMap[IndMax_s]) /(RealT)(255-IndMax_s);
		for (IndexC ind_s =0; ind_s<=255-IndMax_s ;ind_s++)
			{	//ColorMap[ind_s+IndMax_s] =(ByteT)(colorStep*ind_s)+ColorMap[IndMax_s];
				ColorMap[ind_s+IndMax_s] =ColorMap[IndMax_s];
				//cout<<"ColorMap["<<ind_s+IndMax_s<<"] = "<<(UIntT)ColorMap[ind_s+IndMax_s]<<endl;
			}
	}
	//cout<<"Last portion interpolated"<<endl;
	m_ColorMap = ColorMap;
	//for(Array1dIterC<ByteT>it (ColorMap); it; it++)	cout<<(UIntT) (*it)<<' ';
	//cout<<endl;
	return 1;
}

bool ColorNormHistShapeC::GetColorMap(const ImageC<ByteRGBValueC> &imSrc, const ImageC<ByteRGBValueC> &imDes)
{
	ImageC<ByteT> imR_s(imSrc.Frame());
    	ImageC<ByteT> imG_s(imSrc.Frame());
    	ImageC<ByteT> imB_s(imSrc.Frame());
    	for( Array2dIter2C<ByteRGBValueC, ByteT> itR(imSrc, imR_s); itR;itR++)
    		itR.Data2() = itR.Data1().Red(); 
        for( Array2dIter2C<ByteRGBValueC, ByteT> itG(imSrc, imG_s); itG;itG++)
    		 itG.Data2() = itG.Data1().Green(); 
        for( Array2dIter2C<ByteRGBValueC, ByteT> itB(imSrc, imB_s); itB;itB++)
    		 itB.Data2() = itB.Data1().Blue(); 
    	
	ImageC<ByteT> imR_d(imDes.Frame());
    	ImageC<ByteT> imG_d(imDes.Frame());
    	ImageC<ByteT> imB_d(imDes.Frame());
    	for( Array2dIter2C<ByteRGBValueC, ByteT> itR(imDes, imR_d); itR;itR++)
    		itR.Data2() = itR.Data1().Red(); 
        for( Array2dIter2C<ByteRGBValueC, ByteT> itG(imDes, imG_d); itG;itG++)
    		 itG.Data2() = itG.Data1().Green(); 
        for( Array2dIter2C<ByteRGBValueC, ByteT> itB(imDes, imB_d); itB;itB++)
    		 itB.Data2() = itB.Data1().Blue(); 

	GetColorMap(imR_s, imR_d, m_colorMapR);
	GetColorMap(imG_s, imG_d, m_colorMapG);
	GetColorMap(imB_s, imB_d, m_colorMapB);
    	
	return 1;
};
bool ColorNormHistShapeC::GetColorMap(const ImageC<ByteRGBValueC>&imSrc,const ImageRectangleC &imRecSrc, const ImageC<ByteRGBValueC> &imDes, const ImageRectangleC &imRecDes)
{	ImageC<ByteRGBValueC> imCropSrc(imSrc, imRecSrc);
	ImageC<ByteRGBValueC> imRegionSrc(imCropSrc.Frame());
	for(	Array2dIter2C<ByteRGBValueC, ByteRGBValueC> it(imCropSrc, imRegionSrc); it;it++)
	{	it.Data2() = it.Data1(); }
	ImageC<ByteRGBValueC> imCropDes(imDes, imRecDes);
	ImageC<ByteRGBValueC> imRegionDes(imCropDes.Frame());
	for(	Array2dIter2C<ByteRGBValueC, ByteRGBValueC> it(imCropDes, imRegionDes); it;it++)
	{	it.Data2() = it.Data1(); }
	
	GetColorMap(imRegionSrc,imRegionDes);
	
	return 1;
	
};
ImageC<ByteRGBValueC> ColorNormHistShapeC::Apply(const ImageC<ByteRGBValueC> &im)
{	ImageC<ByteRGBValueC> imOut(im.Frame()) 	;
	for( Array2dIter2C<ByteRGBValueC, ByteRGBValueC> it(im, imOut); it;it++)
	{	it.Data2().Red() = m_colorMapR[(IndexC)it.Data1().Red()]; 
    		it.Data2().Green() = m_colorMapG[(IndexC)it.Data1().Green()]; 
    		it.Data2().Blue() = m_colorMapB[(IndexC)it.Data1().Blue()]; 
	}
      	return imOut;
   	
};
