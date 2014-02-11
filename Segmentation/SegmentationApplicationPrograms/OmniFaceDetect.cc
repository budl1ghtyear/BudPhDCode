#include "Ravl/Option.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/String.hh"
#include "Ravl/Assert.hh"
#include "Ravl/IO.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ByteRGBValue.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/Image/SubSample.hh"

#include "Omni/ColossusException.hh"
#include "Omni/FaceScan.hh"
#include "Omni/DetectFace.hh"
#include "Omni/FaceDetectionModel.hh"
#include "Omni/DetectedFace.hh"
#include "Omni/RawFace.hh"
#include "Omni/DisplayFace.hh"

using namespace RavlN;
using namespace OmniN;

int main(int argc, char** argv)
{	OptionC opt(argc, argv);
	FilenameC ImageFile = opt.String("i","","Input Query Image");
	FilenameC  faceDetectionModelHQ = opt.String("h","/vol/vssp/localsoft/External/FaceDetect/models/HighQualityFD.abs","High Quality AAM Data"); //
	FilenameC  faceDetectionModelLQ = opt.String("l","/vol/vssp/localsoft/External/FaceDetect/models/LowQualityFaceScan.abs","Low Quality AAM Data"); //
	opt.Check();

	FaceScanC detectAAMFace;           // face locator that can provide face shape from AAM
	DetectFaceC detectFace;            // eezy-use face locator
	DetectedFaceC detectedFace;
	Point2dC chin;
	
  	// set up face and eye locators
    FaceDetectionModelC faceDetectionModel;
    if(!Load(faceDetectionModelHQ, faceDetectionModel))
    	RavlIssueError("Trouble loading eezy-use face detection model\n");
    detectFace = DetectFaceC(faceDetectionModel); 
    if(!Load(faceDetectionModelLQ, detectAAMFace))
    	RavlIssueError("Trouble loading AAM face detection model\n");
      
    ImageC<ByteRGBValueC> fullSizeImg;
    RavlAlwaysAssertMsg(Load(ImageFile, fullSizeImg), "Could not load image: " + ImageFile);	
    UIntT factor(1);
    DetectedFaceC tmp = detectFace.Apply(fullSizeImg);
    
    
    tmp.LeftEye()  *= factor;
    tmp.RightEye() *= factor;
    detectedFace = DetectedFaceC(fullSizeImg, tmp.LeftEye(), tmp.RightEye(), tmp.Quality());
    chin = Point2dC(RavlConstN::nanReal,RavlConstN::nanReal);
    
    Tuple2C<Point2dC,Point2dC> eyes(Point2dC(detectedFace.LeftEye()[0], detectedFace.LeftEye()[1]),Point2dC(detectedFace.RightEye()[0],detectedFace.RightEye()[1])); 
    detectFace.Apply(fullSizeImg, detectedFace);
    //SArray1dC<Point2dC> face_shape = detectAAMFace.GetShape();
    //detectFace.GetFaceShape(face_shape);
    RawFaceC face(fullSizeImg, eyes.Data1(),eyes.Data2());
    DisplayFace(face, "Detected Face");/*
    ByteRGBValueC eyeloc(255,255,255);
    DrawCross(img,eyeloc,eyes.Data1(),5);
    DrawCross(img,eyeloc,eyes.Data2(),5);
    RavlN::Save(("@X: Detected Face and Eyes"), img);
    cerr << "Eyes detected at " << detectedFace.LeftEye() << "; " << detectedFace.RightEye() << endl;*/
    return 1;           
}
