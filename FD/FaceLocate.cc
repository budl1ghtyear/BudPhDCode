#include "Ravl/Image/FaceLocate.hh"
#include "Ravl/StdConst.hh"

#include "Ravl/IO.hh"

#include "Omni/ColossusException.hh"
#include "Omni/FaceDetectionModel.hh"


using namespace RavlImageN;
using namespace OmniN;

namespace RavlN{


  FaceLocateC::FaceLocateC(bool UseLQ)
    : useLQ(UseLQ)
  {
    StringC modelPath = OMNIHOME"/share/FaceDetect/models/";
    if (useLQ) {  // Load slower AAM-base model, for lower quality images
      // Also has possibility of detecting other face parts
      if(!Load(modelPath + "LowQualityFaceScan.abs", detectAAMFace))
        RavlIssueError("Trouble loading LQ face detection model\n");
    }
    else{ // Load faster model, for higher quality images
      FaceDetectionModelC faceDetectionModel;
      if(!Load(modelPath + "HighQualityFD.abs", faceDetectionModel))
        RavlIssueError("Trouble loading HQ face detection model\n");
      detectFace = DetectFaceC(faceDetectionModel); 
    }
  }
  
  void FaceLocateC::Apply(const ImageC<ByteRGBValueC>& img)
  {
    if (useLQ) { // using model for LQ images
      EyePairsListC eyeList;
      DListC<FaceWithEyesRegionsC> faceAndEyesList;
      ImageC<RealT> grayImage(img.Frame());
      for (Array2dIter2C<ByteRGBValueC,RealT> p(img,grayImage); p; ++p) p.Data2() = p.Data1().Y();
      detectAAMFace.Apply(grayImage,  eyeList);

      EyePairC eyes = eyeList.First();
      Point2dC chin = detectAAMFace.GetShape()[7];
      face = DetectedFaceC(img, DetectedEyeC(eyes.lEye[0],eyes.lEye[1]), DetectedEyeC(eyes.rEye[0],eyes.rEye[1]), eyes.reliability);
    }
    else { // using model for HQ images
      face = detectFace.Apply(img);
    }  
  }

}




