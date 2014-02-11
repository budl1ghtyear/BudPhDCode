#ifndef ASMROTATIONNORMALISATION_HEADER
#define ASMROTATIONNORMALISATION_HEADER 1

#include "ASMNormalisation.hh"

using namespace RavlImageN;
using namespace RavlN;

  //! userlevel=Develop
  //: Shape model with scale and rotation normalisation.

  class ASMScaleRotationNormalisationBodyC
    : public ASMNormalisationBodyC
  {
  public:
    ASMScaleRotationNormalisationBodyC();
    //: Default constructor.

    ASMScaleRotationNormalisationBodyC(BinIStreamC &is);
    //: Load from bin stream.

    ASMScaleRotationNormalisationBodyC(istream &is);
    //: Load from stream.

    virtual bool Save(BinOStreamC &out) const;
    //: Save to binary stream 'out'.

    virtual bool Save(ostream &out) const;
    //: Save to stream 'out'.

    virtual bool ComputeMean(SampleStreamC<AAMAppearanceC> &sample);
    //: Compute mean control points for the list of appearance provided.

    virtual bool RawParameters(const AAMAppearanceC &inst,VectorC &fixedParams,VectorC &freeParams) const;
    //: Generate raw parameters.
    //!param: inst        - input appearance for which we would like to compute the parameters.
    //!param: fixedParams - output pose parameters (e.g. pose, scale, orientation).
    //!param: freeParams  - output normalised control point coordinates. This vector consists of the concatenation of the X and Y coordinates of all control points in a normalised frame.
    //  The raw parameters are the parameters representing the shape before applying PCA. They consists of the pose parameters, which describe the pose of the model instance in the image, and the shape parameters, which describe its shape.

    virtual void RawProject(const VectorC &fixedParams,const VectorC &freeParams,SArray1dC<Point2dC> &out) const;
    //: Generate control points defining an appearance from the raw parameters.
    //!param: fixedParams - input pose parameters (e.g. pose, scale, orientation).
    //!param: freeParams  - input normalised control point coordinates. This vector consists of the concatenation of the X and Y coordinates of all control points in a normalised frame.
    //!param: out         - ouput control points

    virtual IntT NoFixedParameters() const;
    //: Return number of parameters describing the pose
    //  These parameters include e.g. the position, scale and orientation of the model instance

  protected:

  };

  //! userlevel=Normal
  //: Shape model with scale and rotation normalisation.

  class ASMScaleRotationNormalisationC
    : public ASMNormalisationC
  {
  public:
    ASMScaleRotationNormalisationC()
    {}
    //: Default constructor.

    ASMScaleRotationNormalisationC(bool)
      : ASMNormalisationC(*new ASMScaleRotationNormalisationBodyC())
    {}
    //: Constructor

    ASMScaleRotationNormalisationC(BinIStreamC &is);
    //: Binary stream constructor.

    ASMScaleRotationNormalisationC(istream &is);
    //: Stream constructor.

  protected:
    ASMScaleRotationNormalisationC(ASMNormalisationBodyC &bod)
      : ASMNormalisationC(bod)
    {}
    //: Body constructor.

    ASMScaleRotationNormalisationC(ASMScaleRotationNormalisationBodyC *bod)
      : ASMNormalisationC(bod)
    {}
    //: Body ptr constructor.

    ASMScaleRotationNormalisationBodyC &Body()
    { return static_cast<ASMScaleRotationNormalisationBodyC &>(ASMNormalisationC::Body()); }
    //: Access body.

    const ASMScaleRotationNormalisationBodyC &Body() const
    { return static_cast<const ASMScaleRotationNormalisationBodyC &>(ASMNormalisationC::Body()); }
    //: Access body.

  public:
  };

  inline
  BinOStreamC &operator<<(BinOStreamC &s,const ASMScaleRotationNormalisationC &ap) {
    ap.Save(s);
    return s;
  }
  //: Save shape model with scale and rotation normalisation to binary stream.

  inline
  BinIStreamC &operator>>(BinIStreamC &s,ASMScaleRotationNormalisationC &ap) {
    ap = ASMScaleRotationNormalisationC(s);
    return s;
  }
  //: Read shape model with scale and rotation normalisation from binary stream.


#endif
