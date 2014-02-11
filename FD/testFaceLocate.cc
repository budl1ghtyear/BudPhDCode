
#include <cstdlib>
#include "Ravl/Option.hh"
#include "Ravl/IO.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/DP/SequenceIO.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ByteRGBValue.hh"
#include "Ravl/Image/SubSample.hh"
//#include "Ravl/Image/FaceLocate.hh"
#include "FaceLocate.hh"
#include "Omni/ColossusException.hh"
#include "Omni/FaceScan.hh"
#include "Omni/DetectFace.hh"
#include "Omni/FaceDetectionModel.hh"
#include "Omni/DetectedFace.hh"
#include "Ravl/VirtualConstructor.hh"

using namespace RavlN;
using namespace RavlImageN;
using namespace OmniN;

int main (int argc, char* argv[]) {
  //CreateVirtualConstructorAlias("AAMTScaleRotationShapeModelBodyC","AAMScaleRotationShapeModelBodyC");
  ImageC<ByteRGBValueC> img;
  Load(StringC(getenv("PROJECT_OUT")) + "/share/RAVL/testData/face.jpg", img);

  RavlN::FaceLocateC faceLocate(false);
  faceLocate.Apply(img);
  if ((   faceLocate.Quality() < 0.99)
      || (Abs(faceLocate.RightEye().Col() - 312.4) > 1)
      || (Abs(faceLocate.RightEye().Row() - 279.2) > 1)
      || (Abs(faceLocate.LeftEye() .Col() - 416.6) > 1)
      || (Abs(faceLocate.LeftEye(). Row() - 280.1) > 1)
      ) {
    cout << "Error in HQ face location\n";
    cout << faceLocate.RightEye()<<endl;
    cout << faceLocate.LeftEye()<<endl;
    return __LINE__;
  }

  /*
  RavlN::FaceLocateC faceLocateAAM(true);
  faceLocateAAM.Apply(img);
  if ((   faceLocate.Quality() < 0.99)
      || (Abs(faceLocateAAM.RightEye().Col() - 312.4) > 2)
      || (Abs(faceLocateAAM.RightEye().Row() - 279.2) > 2)
      || (Abs(faceLocateAAM.LeftEye() .Col() - 416.6) > 2)
      || (Abs(faceLocateAAM.LeftEye(). Row() - 280.1) > 2)
      ) {
    cout << "Error in LQ face location\n";
    cout << faceLocateAAM.Quality()<<endl;
    cout << faceLocateAAM.RightEye()<<endl;
    cout << faceLocateAAM.LeftEye()<<endl;
    return __LINE__;
  }
  */
}

