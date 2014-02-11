#ifndef RAVL_FACELOCATE_HH
#define RAVL_FACELOCATE_HH

//! author="Bill Christmas"
//! docentry="Ravl.API.Images.Face_Detection"
//! lib=RavlImageProc

#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ByteRGBValue.hh"
#include "Omni/DetectedFace.hh"
#include "Omni/FaceScan.hh"
#include "Omni/DetectFace.hh"



using namespace RavlImageN;
using namespace OmniN;

namespace RavlN {

  //! userlevel=Normal

  //: Class to generate eye positions from frontal face images.
  // <p>This face locator uses some of the OmniPerception algorithms to generate
  // the positions of the eye socket centres.  It is designed <i>only</i> for images
  // that are more or less frontal.</p>
  // <p>Two algorithms are available:</p><ul>

  // <li> A relatively fast one, intended for use with high quality frontal
  // images.</li>

  // <li> A slower one that should handle lower quality (but still frontal)
  // images.</li></ul>

  // <p>Currently this class is only available for shared library linking
  // options (qm shared etc.)</p>

  class FaceLocateC {
  public:
    
    FaceLocateC(bool UseLQ=false);
    //: Constructor
    // Selects and loads the appropriate face model, as described above
    //!param: UseLQ - uses LQ model if true, HQ model otherwise (the default)
    
    void Apply(const ImageC<ByteRGBValueC>& img);
    //: Computes the eye locations in image "img"
    
    const Point2dC LeftEye()
    { return face.LeftEye(); }
    //: Returns the left eye of the person
    // i.e. the eye to the right in the image
    
    const Point2dC RightEye()
    { return face.RightEye(); }
    //: Returns the right eye of the person
    // i.e. the eye to the left in the image
    
    const RealT Quality() const
    { return face.Quality(); }
    //: Returns a number in range (0..1) indicating confidence in the result
    
  private:

    bool useLQ;
    DetectedFaceC face;
    FaceScanC detectAAMFace;           
    DetectFaceC detectFace;

  };

}

#endif
