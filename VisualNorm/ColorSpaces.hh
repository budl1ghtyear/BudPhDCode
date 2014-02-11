#ifndef COLORSPACES_HEADER
#define COLORSPACES_HEADER 1

#include "Ravl/Vector.hh"
#include "Ravl/Vector3d.hh"
#include "Ravl/Matrix3d.hh"
#include "math.h"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ByteRGBValue.hh"
#include "Ravl/Image/RealRGBValue.hh"

using namespace RavlN;
using namespace RavlImageN;

Vector3dC rgb2xyz(Vector3dC rgb);
Vector3dC xyz2rgb(Vector3dC xyz);
Vector3dC xyz2Lab(Vector3dC xyz);
Vector3dC Lab2xyz(Vector3dC Lab);
Vector3dC Luv2xyz(Vector3dC Luv);
Vector3dC xyz2Luv(Vector3dC xyz);

RealT PerceptColorDiff(Vector3dC Labstd, Vector3dC Labsample, Vector3dC KLCH);

ImageC<RealRGBValueC> imRGB2XYZ(const ImageC<ByteRGBValueC> &imRGB);
ImageC<RealRGBValueC> imXYZ2LAB(const ImageC<RealRGBValueC> &imXYZ);
ImageC<RealRGBValueC> imLAB2XYZ(const ImageC<RealRGBValueC> &imLAB);
ImageC<ByteRGBValueC> imXYZ2RGB(const ImageC<RealRGBValueC> &imXYZ);
ImageC<RealRGBValueC> imLUV2XYZ(const ImageC<RealRGBValueC> &imLUV);
ImageC<RealRGBValueC> imXYZ2LUV(const ImageC<RealRGBValueC> &imXYZ);

#endif
