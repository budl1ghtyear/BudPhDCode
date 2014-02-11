#include "Ravl/Option.hh"
#include "Ravl/IO.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/DP/SequenceIO.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ByteRGBValue.hh"
#include "Ravl/Image/SubSample.hh"

#include "Omni/ColossusException.hh"
#include "Omni/FaceScan.hh"
#include "Omni/DetectFace.hh"
#include "Omni/FaceDetectionModel.hh"
#include "Omni/DetectedFace.hh"

using namespace RavlN;
using namespace RavlImageN;
using namespace OmniN;

int main (int argc, char* argv[]) {
  OptionC opt(argc,argv);
  StringC inSeq = opt.String("i", "", "Input face sequence");
  StringC outFile = opt.String("o", "-", "Output file name");
  bool debug = opt.Boolean("d", "debug");
  bool aam  = opt.Boolean("lq", "Use alternative AAM-based model for lower-quality images");
  opt.Compulsory("i");
  opt.Check();

  StringC modelPath = "/vol/vssp/localsoft/External/FaceDetect/models/";
  if (debug) cerr << "Loading face detection models from files in " << modelPath << endl;
  DetectFaceC detectFace;
  FaceScanC detectAAMFace;           

  if (aam) {  // Load slower AAM-base model, for lower quality images
    // Also has possibility of detecting other face parts
    if(!Load(modelPath + "LowQualityFaceScan.abs", detectAAMFace))
      RavlIssueError("Trouble loading AAM face detection model\n");
  }
  else{ // Load faster model, for higher quality images
    FaceDetectionModelC faceDetectionModel;
    if(!Load(modelPath + "HighQualityFD.abs", faceDetectionModel))
      RavlIssueError("Trouble loading eezy-use face detection model\n");
    detectFace = DetectFaceC(faceDetectionModel); 
  }

  DPIPortC<ImageC<ByteRGBValueC> > in;
  if(!OpenISequence(in, inSeq)) { // Failed to open input file. Report an error...
    exit(-1);
  }
  OStreamC faceOut(outFile);
  faceOut << "filename: " << inSeq << "\n"
          << "frames:\n{\n";

  ImageC<ByteRGBValueC> fullSizeImg;
  ImageC<ByteRGBValueC> img;
  while(in.Get(fullSizeImg)) {

    // subsample if too big
    UIntT factor(1);
    if (fullSizeImg.Rows() > 1000) {
      img = fullSizeImg.Copy();
      while (img.Rows() > 1000) {
        img = SubSample(img,2);
        factor *= 2;
      }
    }
    else img = fullSizeImg;
    if (debug) {
      cerr << "Subsampling factor = " << factor << endl;
    }
    

    DetectedFaceC face;
    if (aam) { // using model for LQ images
      EyePairsListC eyeList;
      DListC<FaceWithEyesRegionsC> faceAndEyesList;
      ImageC<RealT> grayImage(img.Frame());
      for (Array2dIter2C<ByteRGBValueC,RealT> p(img,grayImage); p; ++p) p.Data2() = p.Data1().Y();
      detectAAMFace.MFApply(grayImage, FACE_SCAN_SVM, EYE_SCAN_COMBINED_D, 1, true, false, eyeList, faceAndEyesList);
      EyePairC eyes = eyeList.First();
      Point2dC chin = detectAAMFace.GetShape()[7]*factor;
      face = DetectedFaceC(fullSizeImg, DetectedEyeC(eyes.lEye[0]*factor,eyes.lEye[1]*factor), DetectedEyeC(eyes.rEye[0]*factor,eyes.rEye[1]*factor), 0.1);
    }
    else { // using model for HQ images
      face = detectFace.Apply(img);
      face.LeftEye()  *= factor;
      face.RightEye() *= factor;
		if(!Save("@X:Original Face Image",img)) cerr<<"Could not show saved image"<<endl;
		if(!Save("@X:Detected Face Image()",face.Image())) cerr<<"Could not show saved image"<<endl;
		if(!Save("@X:Detected Face Double Image()",face.DoubleImage())) cerr<<"Could not show saved image"<<endl;
    }
    faceOut << "  frame:\n  {\n    face:\n    {\n"
            << "      score: " << face.Quality() <<"\n"
            << "      points:\n      {\n"
            << "        pt: { "
            << "x: " << face.RightEye().Col() << " y: " << face.RightEye().Row()
            << " score: 1 }\n"
            << "        pt: { "
            << "x: " << face.LeftEye().Col() << " y: " << face.LeftEye().Row()
            << " score: 1 }\n"
            << "      }\n    }\n  }\n";
  }
  faceOut << "}\n";
  
}

