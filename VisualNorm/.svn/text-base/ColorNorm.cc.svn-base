////////////////////////////////////////////////////////////////////
//! file="ColorNorm.cc"
//! author="Xuan Zou"

#include         "ColorNorm.hh"
#include 	"Ravl/Vector3d.hh"
#include 	"Ravl/Array2dIter.hh"
#include 	"Ravl/Array2dIter2.hh"
#include 	"Ravl/LeastSquares.hh"
#include 	"Ravl/Array1dIter.hh"
#include 	"Ravl/SArray1dIter.hh"

using namespace RavlN;
using namespace RavlImageN;

VectorC ColorNormC::IllumEstimate(const ImageC<ByteRGBValueC> &im, bool verbose )
{  RealT size =(RealT)im.Size(); 
    //cout<<"pixel number: "<<size<<endl;
    //get R, G ,B channels of the input image
    ImageC<RealT> imR(im.Frame());
    ImageC<RealT> imG(im.Frame());
    ImageC<RealT> imB(im.Frame());
    for( Array2dIter2C<ByteRGBValueC, RealT> itR(im, imR); itR;itR++)
    { itR.Data2() = itR.Data1().Red(); 
    }
    for( Array2dIter2C<ByteRGBValueC, RealT> itG(im, imG); itG;itG++)
    { itG.Data2() = itG.Data1().Green(); 
    }
    for( Array2dIter2C<ByteRGBValueC, RealT> itB(im, imB); itB;itB++)
    { itB.Data2() = itB.Data1().Blue(); 
    }

    //
    switch(m_type)
    {  case 0: //GREY_WORLD
            {  RealT avgR = imR.Sum()/size;
                RealT avgG = imG.Sum()/size;
		RealT avgB = imB.Sum()/size;
		//cout<<"avgR: "<<avgR<<" avgG: "<<avgG<<" avgB: "<<avgB<<endl;
                
		RealT illumR = avgR;
		RealT illumG = avgG;
		RealT illumB = avgB;
		return VectorC(illumR, illumG, illumB);
             }
             break;
        case 1: //CHANNEL_MAX 
            { RealT maxR = imR[imR.IndexOfMax()] ;
               RealT maxG = imG[imG.IndexOfMax()] ;
               RealT maxB = imB[imB.IndexOfMax()] ;
	      
	        RealT illumR = maxR;
		RealT illumG = maxG;
		RealT illumB = maxB;	
               return VectorC(illumR, illumG, illumB);
             
            }
            break;
        case 2:  //GREY_EDGE
            { // to fill
		return VectorC(255/3, 255/3, 255/3);
            }break;
        default:;
    }
    return VectorC(255/3, 255/3, 255/3);
}
VectorC ColorNormC::IllumEstimate(const ImageC<ByteRGBValueC> &im, const ImageRectangleC &imRec , bool verbose )
{	ImageC<ByteRGBValueC> imCrop(im, imRec);
	ImageC<ByteRGBValueC> imRegion(imCrop.Frame());
	for(	Array2dIter2C<ByteRGBValueC, ByteRGBValueC> it(imCrop, imRegion); it;it++)
	{	it.Data2() = it.Data1(); }
	return (IllumEstimate(imRegion, verbose));
}

ImageC<ByteRGBValueC> ColorNormC::Apply(const ImageC<ByteRGBValueC> &im, bool verbose )
{ 
    ImageC<ByteRGBValueC> imOut(im.Frame());
    RealT size =(RealT)im.Size();
    RealT illumR=1, illumG=1, illumB=1;
    //get R, G ,B channels of the input image
    ImageC<RealT> imR(im.Frame());
    ImageC<RealT> imG(im.Frame());
    ImageC<RealT> imB(im.Frame());
    for( Array2dIter2C<ByteRGBValueC, RealT> itR(im, imR); itR;itR++)
    { itR.Data2() = itR.Data1().Red(); 
    }
    for( Array2dIter2C<ByteRGBValueC, RealT> itG(im, imG); itG;itG++)
    { itG.Data2() = itG.Data1().Green(); 
    }
    for( Array2dIter2C<ByteRGBValueC, RealT> itB(im, imB); itB;itB++)
    { itB.Data2() = itB.Data1().Blue(); 
    }
    //estimate the illumination
    switch(m_type)
    {  case 0: //GREY_WORLD
            {  RealT avgR = imR.Sum()/size;
                RealT avgG = imG.Sum()/size;
		RealT avgB = imB.Sum()/size;
                illumR = avgR;
		illumG = avgG;
		illumB = avgB;
	   }
             break;
        case 1: //CHANNEL_MAX 
            { RealT maxR = imR[imR.IndexOfMax()] ;
               RealT maxG = imG[imG.IndexOfMax()] ;
               RealT maxB = imB[imB.IndexOfMax()] ;
	       illumR = maxR;
		illumG = maxG;
		illumB = maxB;	
               }
            break;
        case 2:  //GREY_EDGE
            { 
            }break;
        default:;
    }
    VectorC illum(illumR,illumG, illumB);
    //cout<<"illum :"<<illum<<endl;
    RealT invIllumR= 1/illumR;
    RealT invIllumG= 1/illumG;
    RealT invIllumB= 1/illumB;
    imR = imR * invIllumR;
    imG = imG * invIllumG;
    imB = imB * invIllumB;
    //scale all channels so that the maxim value of all channels reaches 255
    RealT maxR = imR[imR.IndexOfMax()] ;
    RealT maxG = imG[imG.IndexOfMax()] ;
    RealT maxB = imB[imB.IndexOfMax()] ;
    RealT maxC=VectorC(maxR, maxG, maxB).MaxValue();
    RealT factor = 255/maxC;
    //cout<<"Max R, G, B:" <<maxR<<'  '<<maxG<<'  '<<maxB<<'  '<<endl;
    //cout<<"factor: "<<'  '<<factor<<endl;
    imR = imR *factor;
    imG =imG * factor;
    imB= imB *factor;
    
    //assign values to R, G ,B channels of the output image
    for( Array2dIter2C<ByteRGBValueC, RealT> itR(imOut, imR); itR;itR++)
    { itR.Data1().Red() =(ByteT) itR.Data2() ; 
    }
    for( Array2dIter2C<ByteRGBValueC, RealT> itG(imOut, imG); itG;itG++)
    { itG.Data1().Green() = (ByteT)itG.Data2() ; 
    }
    for( Array2dIter2C<ByteRGBValueC, RealT> itB(imOut, imR); itB;itB++)
    {  itB.Data1().Blue()= (ByteT)itB.Data2() ; 
    }

    return imOut;
}

ImageC<ByteRGBValueC> ColorNormC::Transform(const ImageC<ByteRGBValueC> &im, const VectorC illumSrc, const VectorC illumDes, bool verbose  )
{	
    if(verbose) cout<<" destination illum: "<< illumDes[0]<<' '<<illumDes[1]<<' '<<illumDes[2]<<endl;
    ImageC<ByteRGBValueC> imOut(im.Frame());
    
    //get R, G ,B channels of the input image
    ImageC<RealT> imR(im.Frame());
    ImageC<RealT> imG(im.Frame());
    ImageC<RealT> imB(im.Frame());
    for( Array2dIter2C<ByteRGBValueC, RealT> itR(im, imR); itR;itR++)
    { itR.Data2() = itR.Data1().Red(); 
    }
    for( Array2dIter2C<ByteRGBValueC, RealT> itG(im, imG); itG;itG++)
    { itG.Data2() = itG.Data1().Green(); 
    }
    for( Array2dIter2C<ByteRGBValueC, RealT> itB(im, imB); itB;itB++)
    { itB.Data2() = itB.Data1().Blue(); 
    }
    
    if(illumSrc[0] <0.001 || illumSrc[1] <0.001 || illumSrc[2]<0.001) {cerr<<"zero element in source illum vec"<<endl; exit(0);}	

    RealT invIllumR= 1/illumSrc[0];
    RealT invIllumG= 1/illumSrc[1];
    RealT invIllumB= 1/illumSrc[2];
    imR = imR * invIllumR * illumDes[0];
    imG = imG * invIllumG * illumDes[1];
    imB = imB * invIllumB * illumDes[2];
    //assign values to R, G ,B channels of the output image
    for( Array2dIter2C<ByteRGBValueC, RealT> itR(imOut, imR); itR;itR++)
    { if( itR.Data2()>255)  itR.Data1().Red() = 255 ; 
      else if( itR.Data2()<0) itR.Data1().Red() =0; 
	      else   itR.Data1().Red() =(ByteT)itR.Data2();
    }
    for( Array2dIter2C<ByteRGBValueC, RealT> itG(imOut, imG); itG;itG++)
    { if( itG.Data2()>255)  itG.Data1().Green() = 255 ; 
      else if( itG.Data2()<0) itG.Data1().Green() =0; 
	      else   itG.Data1().Green() =(ByteT)itG.Data2();
    }
    for( Array2dIter2C<ByteRGBValueC, RealT> itB(imOut, imB); itB;itB++)
    {  if( itB.Data2()>255)  itB.Data1().Blue() = 255 ; 
      else if( itB.Data2()<0) itB.Data1().Blue() =0; 
	      else   itB.Data1().Blue() =(ByteT)itB.Data2();
    }
    // cout<<"ImageOut[1]: " <<imOut[1]<<endl;
    return imOut;
}


ImageC<ByteRGBValueC> ColorNormC::Transform(const ImageC<ByteRGBValueC> &im, const VectorC illumDes, bool verbose  )
{  if(verbose) cout<<" destination illum: "<< illumDes[0]<<' '<<illumDes[1]<<' '<<illumDes[2]<<endl;
    ImageC<ByteRGBValueC> imOut(im.Frame());
    RealT size =(RealT)im.Size();
    RealT illumR=1, illumG=1, illumB=1;
    //get R, G ,B channels of the input image
    ImageC<RealT> imR(im.Frame());
    ImageC<RealT> imG(im.Frame());
    ImageC<RealT> imB(im.Frame());
    for( Array2dIter2C<ByteRGBValueC, RealT> itR(im, imR); itR;itR++)
    { itR.Data2() = itR.Data1().Red(); 
    }
    for( Array2dIter2C<ByteRGBValueC, RealT> itG(im, imG); itG;itG++)
    { itG.Data2() = itG.Data1().Green(); 
    }
    for( Array2dIter2C<ByteRGBValueC, RealT> itB(im, imB); itB;itB++)
    { itB.Data2() = itB.Data1().Blue(); 
    }
    //estimate the illumination
    switch(m_type)
    {  case 0: //GREY_WORLD
            {  RealT avgR = imR.Sum()/size;
                RealT avgG = imG.Sum()/size;
		RealT avgB = imB.Sum()/size;
            	illumR = avgR;
		illumG = avgG;
		illumB = avgB;
	     }
             break;
        case 1: //CHANNEL_MAX 
            { RealT maxR = imR[imR.IndexOfMax()] ;
               RealT maxG = imG[imG.IndexOfMax()] ;
               RealT maxB = imB[imB.IndexOfMax()] ;
	  	illumR = maxR;
		illumG = maxG;
		illumB = maxB;	
             }
            break;
        case 2:  //GREY_EDGE
            { 
            }break;
        default:;
    }
    VectorC illum(illumR,illumG, illumB);
    if (verbose) cout<<"illum of current image:"<<illum<<endl;

    RealT invIllumR= 1/illumR;
    RealT invIllumG= 1/illumG;
    RealT invIllumB= 1/illumB;
    imR = imR * invIllumR * illumDes[0];
    imG = imG * invIllumG * illumDes[1];
    imB = imB * invIllumB * illumDes[2];
   
    //assign values to R, G ,B channels of the output image
    for( Array2dIter2C<ByteRGBValueC, RealT> itR(imOut, imR); itR;itR++)
    { if( itR.Data2()>255)  itR.Data1().Red() = 255 ; 
      else if( itR.Data2()<0) itR.Data1().Red() =0; 
	      else   itR.Data1().Red() =(ByteT)itR.Data2();
    }
    for( Array2dIter2C<ByteRGBValueC, RealT> itG(imOut, imG); itG;itG++)
    { if( itG.Data2()>255)  itG.Data1().Green() = 255 ; 
      else if( itG.Data2()<0) itG.Data1().Green() =0; 
	      else   itG.Data1().Green() =(ByteT)itG.Data2();
    }
    for( Array2dIter2C<ByteRGBValueC, RealT> itB(imOut, imB); itB;itB++)
    {  if( itB.Data2()>255)  itB.Data1().Blue() = 255 ; 
      else if( itB.Data2()<0) itB.Data1().Blue() =0; 
	      else   itB.Data1().Blue() =(ByteT)itB.Data2();
    }
    // cout<<"ImageOut[1]: " <<imOut[1]<<endl;
    return imOut;
}

bool ColorNormAffineC::LearnTransformFromImgs(const ImageC<ByteRGBValueC>& imSrc, const ImageC<ByteRGBValueC>& imDes,bool verbose )
{  if ( imSrc.Cols()!= imDes.Cols() || imSrc.Rows() !=imDes.Rows() ) 
    { cerr<<"Dimensions for images do not match! "<<endl; return 0; }
    ImageC<RealRGBValueC> imSrc_R(imSrc.Frame());
    ImageC<RealRGBValueC> imDes_R(imSrc.Frame());
    for(Array2dIter2C<ByteRGBValueC, RealRGBValueC> itSrc(imSrc,imSrc_R);	itSrc;itSrc++)
	{	itSrc.Data2().Red() = itSrc.Data1().Red();
		itSrc.Data2().Green() = itSrc.Data1().Green();
		itSrc.Data2().Blue() = itSrc.Data1().Blue();
	}
    for(Array2dIter2C<ByteRGBValueC, RealRGBValueC> itDes(imDes,imDes_R); itDes;itDes++)
	{	itDes.Data2().Red() = itDes.Data1().Red();
		itDes.Data2().Green() = itDes.Data1().Green();
		itDes.Data2().Blue() = itDes.Data1().Blue();
	}	
    Affine3dC affine = LearnTransform(imSrc_R, imDes_R,verbose);
    SetTransform(affine);
    return 1;
}
bool ColorNormAffineC::LearnTransformFromImgs(const ImageC<RealRGBValueC>& imSrc, const ImageC<RealRGBValueC>& imDes,bool verbose )
{  if ( imSrc.Cols()!= imDes.Cols() || imSrc.Rows() !=imDes.Rows() ) 
    { cerr<<"Dimensions for images do not match! "<<endl; return 0; }
    Affine3dC affine = LearnTransform(imSrc, imDes,verbose);
    SetTransform(affine);
    return 1;
}

Affine3dC ColorNormAffineC::LearnTransform(const ImageC<RealRGBValueC>& imSrc, const ImageC<RealRGBValueC>& imDes, bool verbose )
{  if ( imSrc.Cols()!= imDes.Cols() || imSrc.Rows() !=imDes.Rows() ) 
    { cerr<<"Dimensions for images do not match! "<<endl; exit(0);}

    VectorC vecTrans(12);
    vecTrans.Fill(0);
    
    IntT rows = imSrc.Cols();
    IntT cols = imSrc.Rows();
    MatrixC matSrc( 3*rows*cols, 12);
    matSrc.Fill(0);
    UIntT index = 0;
    for(Array2dIterC<RealRGBValueC> it(imSrc) ; it;it++)
    {	matSrc[index][0] = (*it).Red(); matSrc[index][1]=(*it).Green();  matSrc[index][2]=(*it).Blue();  matSrc[index][3]=1;
        matSrc[index+1][4] = (*it).Red(); matSrc[index+1][5]=(*it).Green();  matSrc[index+1][6]=(*it).Blue();  matSrc[index+1][7]=1;
	matSrc[index+2][8] = (*it).Red(); matSrc[index+2][9]=(*it).Green();  matSrc[index+2][10]=(*it).Blue();  matSrc[index+2][11]=1;
        index=index+3;
    }

    VectorC vecDes(3*rows*cols);
    vecDes.Fill(0);
    index = 0;
    for(Array2dIterC<RealRGBValueC> it(imDes) ; it;it++)
    {	vecDes[index] = (*it).Red(); vecDes[index+1]=(*it).Green();   vecDes[index+2]=(*it).Blue();
        index=index+3;
    }
    
    //LSE solution for matSrc*vecTrans = vecDes
    RealT residual;
    vecTrans = LeastSquaresQR(matSrc, vecDes, residual);
    if(verbose) cout<<"-- Residual :"<<residual<<endl;
    Affine3dC affine(Matrix3dC(vecTrans[0],vecTrans[1],vecTrans[2],vecTrans[4],vecTrans[5],vecTrans[6],vecTrans[8],vecTrans[9],vecTrans[10]) ,
Vector3dC(vecTrans[3], vecTrans[7], vecTrans[11]));
    if(verbose) cout << affine;
    return affine;

}

ImageC<ByteRGBValueC> ColorNormAffineC::Apply(const ImageC<ByteRGBValueC>& im, bool verbose )
{ ImageC<ByteRGBValueC> imOut (im.Frame());
   for(Array2dIter2C<ByteRGBValueC, ByteRGBValueC> it(im, imOut) ; it; it++)
   { Vector3dC vec ; vec[0] = it.Data1().Red(); vec[1] =it.Data1().Green(); vec[2] = it.Data1().Blue(); 
      Vector3dC vecOut = m_affine *vec;
      for( UIntT id=0; id<3;id++)
      { if (vecOut[id] >255)  vecOut[id] =255; else if (vecOut[id]<0) vecOut[id]=0;}
      it.Data2().Red()= (ByteT)vecOut[0]; it.Data2().Green() =(ByteT)vecOut[1]; it.Data2().Blue() =(ByteT)vecOut[2];
   }
   return imOut;
} 
ImageC<RealRGBValueC> ColorNormAffineC::Apply(const ImageC<RealRGBValueC>& im, bool verbose )
{ ImageC<RealRGBValueC> imOut (im.Frame());
   for(Array2dIter2C<RealRGBValueC, RealRGBValueC> it(im, imOut) ; it; it++)
   { Vector3dC vec ; vec[0] = it.Data1().Red(); vec[1] =it.Data1().Green(); vec[2] = it.Data1().Blue(); 
      Vector3dC vecOut = m_affine *vec;
      for( UIntT id=0; id<3;id++)
      { if (vecOut[id] >255)  vecOut[id] =255; else if (vecOut[id]<0) vecOut[id]=0;}
      it.Data2().Red()= vecOut[0]; it.Data2().Green() =vecOut[1]; it.Data2().Blue() =vecOut[2];
   }
   return imOut;
} 

