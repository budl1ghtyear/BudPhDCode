////////////////////////////////////////////////////////////////////
//! file="ColorSpaces.cc"
//! author="Xuan Zou"
//reference: http://www.brucelindbloom.com/
// K. Plataniotis and A. Venetsanopoulos, Color Image Processing and Applications,Springer, 2000
// Danny Pascale, A review of RGB color spaces
#ifndef PI 
#define PI 3.14159
#endif

#include "ColorSpaces.hh"
#include "Ravl/Array2dIter2.hh"


//D65, CIE1964, x y z,   mat*(255,255,255)
Vector3dC xyzD65(0.9505, 1.0000, 1.0890);//xyzD65(0.31382, 0.33100, 0.35518);
Vector3dC xyzD65Inv(1.0521, 1, 0.9183); //xyzD65Inv(1/0.31382, 1/0.33100, 1/0.35518);


// Convert from RGB (0-255) to xyz 
Vector3dC rgb2xyz(Vector3dC rgb)
{ 	// sRGB,  Reference White: D65 
	Matrix3dC mat( 0.412424,  0.357579,  0.180464,  0.212656,  0.715158,  0.0721856,  0.0193324,  0.119193, 0.950444 );
	return mat*(rgb/255);
}
// convert normalised xyz( x+y+z=1) to RGB
Vector3dC xyz2rgb(Vector3dC xyz)
{ 	// sRGB,  Reference White: D65 
	Matrix3dC mat( 3.24071,  -1.53726 ,  -0.498571,   -0.969258 ,  1.87599, 0.0415557, 0.0556352,  -0.203996,  1.05707 );
	return mat*xyz*255;
}


Vector3dC xyz2Lab(Vector3dC xyz)
{	RealT yyn = xyz[1]*xyzD65Inv[1];
        RealT temp= pow(yyn, 0.333);
	RealT L; 
        if(yyn>0.008856) L = 116*temp -16 ;
	else L=903.3*yyn;
	RealT a = 500*(pow(xyz[0]*xyzD65Inv[0], 0.333) - temp);
	RealT b = 200*(temp - pow(xyz[2]*xyzD65Inv[2], 0.333));
	return Vector3dC(L,a,b);
} 

Vector3dC Lab2xyz(Vector3dC Lab)
{	RealT x,y,z; 
        if( Lab[0] >8.00) 
	{ x = xyzD65[0] * pow( (Lab[0]+16)/116 + Lab[1]/500, 3);
	   y = xyzD65[1] * pow( (Lab[0]+16)/116 , 3);
	   z = xyzD65[2] * pow( (Lab[0]+16)/116 - Lab[2]/200, 3);	
	}
         else
	{ x = xyzD65[0] * pow( pow(Lab[0]/903.3 , 0.333) + Lab[1]/500, 3);
	   y = xyzD65[1] * (Lab[0]/903.3);
	   z = xyzD65[2] * pow( pow(Lab[0]/903.3 , 0.333) - Lab[2]/200, 3);	
	}
	return Vector3dC(x,y,z);
}
//ref: K. Plataniotis and A. Venetsanopoulos, Color Image Processing and Applications,Springer, 2000
Vector3dC xyz2Luv(Vector3dC xyz)
{	if(xyz[0]<0.0001 &&xyz[1]<0.0001 && xyz[2]<0.0001)
 	  { xyz[0]=0.0001; xyz[1]=0.0001; xyz[2]=0.0001;}
	RealT yyn = xyz[1]*xyzD65Inv[1];
        RealT temp= pow(yyn, 0.333);
	RealT L; 
        if(yyn>0.008856) L = 116*temp -16 ;
	else L=903.3*yyn;
	RealT tempxyz = 	xyz[0]+15*xyz[1]+3*xyz[2];
	RealT tempxyzn = xyzD65[0]+15*xyzD65[1]+3*xyzD65[2];
	RealT uprime = 4*xyz[0]/tempxyz;
	RealT unprime = 4*xyzD65[0]/tempxyzn;
	RealT vprime = 9*xyz[1]/tempxyz;
	RealT vnprime = 9*xyzD65[1]/tempxyzn;
	RealT u = 13*L*(uprime - unprime);
	RealT v = 13*L*(vprime - vnprime);
	return Vector3dC(L,u,v);
} 
Vector3dC Luv2xyz(Vector3dC Luv)
{	RealT x,y,z; 
        if( Luv[0] >8.00) 
		y = xyzD65[1]*pow((Luv[0]+16) /116, 3);
	else
		y = xyzD65[1]*Luv[0]/903.3;
	RealT  tempxyzn = xyzD65[0]+15*xyzD65[1]+3*xyzD65[2];
	RealT unprime = 4*xyzD65[0]/tempxyzn;
	RealT vnprime = 9*xyzD65[1]/tempxyzn;
	RealT uprime = Luv[1]/(13*Luv[0])+ unprime;
	RealT vprime = Luv[2]/(13*Luv[0])+ vnprime;
	x = 0.25*(uprime *y*(9.0-15.0*vprime)/vprime + 15.0*uprime*y);
	z = ((9.0-15.0*vprime)*y/vprime - x)/3;
	/*RealT a = (52*Luv[0]/(Luv[1]+13*Luv[0]*unprime)-1)/3;
	RealT b = -5*y;
	RealT c= -1/3;
	RealT d =y*(39*Luv[0] /(Luv[2]+13*Luv[0]*vnprime)-5);
	x = (d-b)/(a-c);
	z = a*x+b;
	*/
	return Vector3dC(x,y,z);
}

// Compute the CIEDE2000 color-difference between the sample between a reference
//with CIELab coordinates Labsample and a standard with CIELab coordinates 
// Labstd
// The optional argument KLCH is a 1x3 vector containing the
//the value of the parametric weighting factors kL, kC, and kH
// these default to 1. 

// Based on the article:
// "The CIEDE2000 Color-Difference Formula: Implementation Notes, 
// Supplementary Test Data, and Mathematical Observations,", G. Sharma, 
// W. Wu, E. N. Dalal, submitted to Color Research and Application, 
// January 2004.
// available at http://www.ece.rochester.edu/~/gsharma/ciede2000/

RealT PerceptColorDiff(Vector3dC Labstd, Vector3dC Labsample, Vector3dC KLCH)
{	RealT kl =KLCH[0], kc=KLCH[1], kh =KLCH[2];
	RealT Lstd = Labstd[0];
        RealT astd = Labstd[1];
	RealT bstd = Labstd[2];
        RealT Cabstd = sqrt(astd*astd+bstd*bstd);
        RealT Lsample = Labsample[0];
	RealT asample = Labsample[1];
	RealT bsample = Labsample[2];
        RealT Cabsample = sqrt(asample*asample+bsample*bsample);
 	RealT Cabarithmean = (Cabstd + Cabsample)/2;
	RealT G = 0.5* ( 1 - sqrt( pow(Cabarithmean,7)/(pow(Cabarithmean,7) + pow(25,7))));
	RealT apstd = (1+G)*astd; // aprime in paper
	RealT apsample = (1+G)*asample; // aprime in paper
	RealT Cpsample = sqrt( apsample*apsample+bsample*bsample);
	RealT Cpstd = sqrt(apstd*apstd+bstd*bstd);

	// Compute product of chromas and locations at which it is zero for use later
	RealT Cpprod = Cpsample*Cpstd;
	//zcidx = find(Cpprod == 0);

        // Ensure hue is between 0 and 2PI
	RealT hpstd, hpsample;	
	if((fabs(apstd)+fabs(bstd))<0.00001) hpstd =0;
	else 
	{ hpstd = atan2(bstd,apstd);
           if(hpstd<0) hpstd = hpstd+2*PI; // rollover ones that come -ve
	}
	if((fabs(apsample)+fabs(bsample))<0.00001) hpsample =0;
	else 
	{ hpsample = atan2(bsample,apsample);
           if(hpsample<0) hpsample = hpsample+2*PI; // rollover ones that come -ve
	}
	
	RealT dL = (Lsample-Lstd);
	RealT dC = (Cpsample-Cpstd);
       // Computation of hue difference
	RealT dhp = (hpsample-hpstd);
	if(fabs(Cpsample*Cpstd) <0.000001) dhp = 0;
        else if (dhp>PI)  dhp = dhp - 2*PI;
	else if (dhp< -PI)  dhp = dhp + 2*PI;

	// set chroma difference to zero if the product of chromas is zero

	// Note that the defining equations actually need
	// signed Hue and chroma differences which is different
	// from prior color difference formulae

	RealT dH = 2*sqrt(Cpprod)*sin(dhp/2);
	//dH2 = 4*Cpprod.*(sin(dhp/2)).^2;

	// weighting functions
	RealT Lp = (Lsample+Lstd)/2;
	RealT Cp = (Cpstd+Cpsample)/2;
	// Average Hue Computation
	// This is equivalent to that in the paper but simpler programmatically.
	//Note average hue is computed in radians and converted to degrees only 
	// where needed
	RealT hp = (hpstd+hpsample)/2;
        if(fabs(Cpsample*Cpstd) <0.00001) hp = hpsample+hpstd;
        // abs hue diff exceeds 180 degrees 
	else if(fabs(hpstd-hpsample)  > PI)
	 	{  if(hpstd+hpsample > 2*PI ) hp = hp - PI; 
	             else  	hp = hp+ PI;
		}
        // Check if one of the chroma values is zero, in which case set 
	// mean hue to the sum which is equivalent to other value
	
	RealT Lpm502 = (Lp-50)*(Lp-50);
	RealT Sl = 1 + 0.015*Lpm502/sqrt(20+Lpm502);  
	RealT Sc = 1+0.045*Cp;
	RealT T = 1 - 0.17*cos(hp - PI/6 ) + 0.24*cos(2*hp) + 0.32*cos(3*hp+PI/30)  -0.20*cos(4*hp-63*PI/180);
	RealT Sh = 1 + 0.015*Cp*T;
	RealT delthetarad = 30*exp(- pow( (180/PI*hp-275)/25, 2));
	RealT Rc =  2*sqrt(pow(Cp,7)/(pow(Cp,7) +pow( 25,7)));
	RealT RT =  - sin(2*delthetarad)*Rc;

	RealT klSl = kl*Sl;
	RealT kcSc = kc*Sc;
	RealT khSh = kh*Sh;
        // The CIE 00 color difference
	RealT de00 = sqrt( pow(dL/(klSl),2) + pow(dC/(kcSc),2) + pow(dH/(khSh),2) + RT*dC/(kcSc)*dH/(khSh) );

     return de00;
}

ImageC<RealRGBValueC> imRGB2XYZ(const ImageC<ByteRGBValueC> &imRGB)
{	ImageC<RealRGBValueC> imXYZ(imRGB.Frame());
	for(Array2dIter2C<RealRGBValueC, ByteRGBValueC> it(imXYZ, imRGB); it; it++)
	{	Vector3dC rgb(it.Data2().Red(), it.Data2().Green(), it.Data2().Blue());
		Vector3dC xyz = rgb2xyz(rgb);
		it.Data1().Red() = xyz[0]; it.Data1().Green() = xyz[1]; it.Data1().Blue() = xyz[2];
	}
	return imXYZ;
}

ImageC<RealRGBValueC> imXYZ2LAB(const ImageC<RealRGBValueC> &imXYZ)
{	ImageC<RealRGBValueC> imLAB(imXYZ.Frame());
	for(Array2dIter2C<RealRGBValueC, RealRGBValueC> it(imLAB, imXYZ); it; it++)
	{	Vector3dC xyz(it.Data2().Red(), it.Data2().Green(), it.Data2().Blue());
		Vector3dC lab = xyz2Lab(xyz);
		it.Data1().Red() = lab[0]; it.Data1().Green() = lab[1]; it.Data1().Blue() = lab[2];
	}
	return imLAB;
}

ImageC<RealRGBValueC> imLAB2XYZ(const ImageC<RealRGBValueC> &imLAB)
{	ImageC<RealRGBValueC> imXYZ(imLAB.Frame());
	for(Array2dIter2C<RealRGBValueC, RealRGBValueC> it(imXYZ, imLAB); it; it++)
	{	Vector3dC lab(it.Data2().Red(), it.Data2().Green(), it.Data2().Blue());
		Vector3dC xyz = Lab2xyz(lab);
		it.Data1().Red() = xyz[0]; it.Data1().Green() = xyz[1]; it.Data1().Blue() = xyz[2];
	}
	return imXYZ;
}

ImageC<ByteRGBValueC> imXYZ2RGB(const ImageC<RealRGBValueC> &imXYZ)
{	ImageC<ByteRGBValueC> imRGB(imXYZ.Frame());
	for(Array2dIter2C<ByteRGBValueC, RealRGBValueC> it(imRGB, imXYZ); it; it++)
	{	Vector3dC xyz(it.Data2().Red(), it.Data2().Green(), it.Data2().Blue());
		Vector3dC rgb = xyz2rgb(xyz);
		for(IntT id = 0;id<3;id++) 
		{	if( rgb[id] <0) rgb[id] =0; else if (rgb[id]>255) rgb[id] = 255;}
		it.Data1().Red() =(ByteT) rgb[0]; it.Data1().Green() =(ByteT) rgb[1]; it.Data1().Blue() = (ByteT) rgb[2];
	}
	return imRGB;
}

ImageC<RealRGBValueC> imXYZ2LUV(const ImageC<RealRGBValueC> &imXYZ)
{	ImageC<RealRGBValueC> imLUV(imXYZ.Frame());
	for(Array2dIter2C<RealRGBValueC, RealRGBValueC> it(imLUV, imXYZ); it; it++)
	{	Vector3dC xyz(it.Data2().Red(), it.Data2().Green(), it.Data2().Blue());
		Vector3dC luv = xyz2Luv(xyz);
		it.Data1().Red() = luv[0]; it.Data1().Green() = luv[1]; it.Data1().Blue() = luv[2];
	}
	return imLUV;
}

ImageC<RealRGBValueC> imLUV2XYZ(const ImageC<RealRGBValueC> &imLUV)
{	ImageC<RealRGBValueC> imXYZ(imLUV.Frame());
	for(Array2dIter2C<RealRGBValueC, RealRGBValueC> it(imXYZ, imLUV); it; it++)
	{	Vector3dC luv(it.Data2().Red(), it.Data2().Green(), it.Data2().Blue());
		Vector3dC xyz = Luv2xyz(luv);
		it.Data1().Red() = xyz[0]; it.Data1().Green() = xyz[1]; it.Data1().Blue() = xyz[2];
	}
	return imXYZ;
}
