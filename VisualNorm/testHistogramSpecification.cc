//////////////////////////////////////////////////////////////////////
//! file="testHistogramSpecification.cc"
//! author="Bud Goswami"
//! this file shows how to use the color transform learnt from a pair of corresponding regions in the source and target image, then apply the transform to the whole image
//! 
#include "Ravl/Point3d.hh"
#include "ColorNorm.hh"
#include "ColorSpaces.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ByteRGBValue.hh"
#include "Ravl/Option.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/OS/Filename.hh" 
#include "Ravl/IO.hh"
#include "Ravl/DP/FileFormatIO.hh"
#include "MouthRegion.hh"

/* color mapping using the transform learnt from a pair of images */
int main(int argc, char **argv)
{ 
   OptionC   opt(argc,argv);
   FilenameC imSrcFile = opt.String("sf","","Source image file for learning the transform" );
   FilenameC imDesFile = opt.String("df","","Target image file for learning the transform " );
   FilenameC imInFile = opt.String("in","","input image file to process" );
   FilenameC imOutFile = opt.String("o","/tmp/temp.ppm","output image file" );
   bool verbose=opt.Boolean("v", false,"verbose");
   opt.Check();
   
   // load images
   ImageC<ByteRGBValueC> origimSrc, origimDes, origimIn, imOut;
   if(!Load(imSrcFile, imSrc))
     { cerr<<"Failed to load image: "<<imSrcFile<<endl;
        return 0;;
     }
   if(!Load(imDesFile, imDes))
     { cerr<<"Failed to load image: "<<imDesFile<<endl;
        return 0;;
     }
   if(!Load(imInFile, imIn))
     { cerr<<"Failed to load image: "<<imInFile<<endl;
        return 0;;
     }
   FilenameC cascade_name = "/opt/share/opencv/haarcascades/haarcascade_frontalface_alt.xml"; 
   ImageRectangleC imIn_face = GetFaceCoords(const origimIn, const FilenameC &c_name);
   ImageC<ByteRGBBValueC> imIn(origimIn,imIn_face);
   Save("@X:Source", imSrc);
   Save("@X:Target", imDes);
   Save("@X:Input", imIn);
  

   ColorNormHistShapeC colorNorm;
	colorNorm.GetColorMap(imSrc, imDes);
	//cout<<"colormap obtained"<<endl;
	imOut = colorNorm.Apply(imIn);
	
   Save("@X: Transformed input", imOut);
  
   if(!Save(imOutFile, imOut))
     { cerr<<"Failed to save image to: "<<imOutFile<<endl;
        return 0;;
     }
 


   return 1;
}
