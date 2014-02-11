#ifndef ASMAFFINENORMALISATION_HEADER
#define ASMAFFINENORMALISATION_HEADER 1

#include "ASMNormalisation.hh"

using namespace RavlImageN ;
using namespace RavlN;

  //! userlevel=Develop
  //: Statistical shape model with affine normalisation.

  class ASMAffineNormalisationBodyC
    : public ASMNormalisationBodyC
  {
  public:
    ASMAffineNormalisationBodyC();
    //: Default constructor.

    ASMAffineNormalisationBodyC(BinIStreamC &is);
    //: Load from bin stream.

    ASMAffineNormalisationBodyC(istream &is);
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
  //: Statistical shape model with affine normalisation.

  class ASMAffineNormalisationC
    : public ASMNormalisationC
  {
  public:
    ASMAffineNormalisationC()
    {}
    //: Default constructor.

    ASMAffineNormalisationC(bool)
      : ASMNormalisationC(*new ASMAffineNormalisationBodyC())
    {}
    //: Constructor

    ASMAffineNormalisationC(BinIStreamC &is);
    //: Binary stream constructor.

    ASMAffineNormalisationC(istream &is);
    //: Stream constructor.

  protected:
    ASMAffineNormalisationC(ASMNormalisationBodyC &bod)
      : ASMNormalisationC(bod)
    {}
    //: Body constructor.

    ASMAffineNormalisationC(ASMAffineNormalisationBodyC *bod)
      : ASMNormalisationC(bod)
    {}
    //: Body ptr constructor.

    ASMAffineNormalisationBodyC &Body()
    { return static_cast<ASMAffineNormalisationBodyC &>(ASMNormalisationC::Body()); }
    //: Access body.

    const ASMAffineNormalisationBodyC &Body() const
    { return static_cast<const ASMAffineNormalisationBodyC &>(ASMNormalisationC::Body()); }
    //: Access body.

  public:
  };

  inline
  BinOStreamC &operator<<(BinOStreamC &s,const ASMAffineNormalisationC &ap) {
    ap.Save(s);
    return s;
  }
  //: Save shape model with affine normalisation to binary stream.

  inline
  BinIStreamC &operator>>(BinIStreamC &s,ASMAffineNormalisationC &ap) {
    ap = ASMAffineNormalisationC(s);
    return s;
  }
  //: Read shape model with affine normalisation from binary stream.



#endif
