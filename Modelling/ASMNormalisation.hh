#ifndef ASMNormalisation_HEADER
#define ASMNormalisation_HEADER 1

#include "Ravl/Image/AAMAppearance.hh"
#include "Ravl/PatternRec/Function.hh"
#include "Ravl/PatternRec/SampleStream.hh"
#include "Ravl/PatternRec/SampleVector.hh"
using namespace RavlImageN;
using namespace RavlN;
  //! userlevel=Develop
  //: Statistical shape model.

  class ASMNormalisationBodyC
    : public RCBodyVC
  {
  public:
    ASMNormalisationBodyC()
      : m_verbose(false), // Write error messages to standard out.
        nPoints(0) // Number of control points in the model.
    {}
    //: Default constructor.

    ASMNormalisationBodyC(BinIStreamC &is);
    //: Constructor.
    //  Load from bin stream.

    ASMNormalisationBodyC(istream &is);
    //: Constructor.
    //  Load from stream.

    void SetVerbose(bool verboseMode)
    { m_verbose = verboseMode; }
    //: If set to true designer will print error messages to std::cerr

    virtual bool Save(BinOStreamC &out) const;
    //: Save to binary stream 'out'.

    virtual bool Save(ostream &out) const;
    //: Save to stream 'out'.

    virtual bool ComputeMean(const SampleC<AAMAppearanceC> &sample);
    //: Compute mean control points for the list of appearance provided.

    virtual bool ComputeMean(SampleStreamC<AAMAppearanceC> &sample);
    //: Compute mean control points for the list of provided appearances

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

    VectorC Parameters(const AAMAppearanceC &inst) const;
    //: Return a parameter vector representing the appearance 'inst'.

    bool Design(const SampleC<AAMAppearanceC> &sample,RealT variation = 0.95,UIntT maxP=50);
    //: Design a shape model given some data.
    //!param: sample    - list of appearances for training the model.
    //!param: variation - percentage of variation preserved during PCA.
    //!param: maxP      - limit on number of parameters contained in the shape model.

    bool Design(SampleStreamC<AAMAppearanceC> &sample,RealT variation = 0.95,UIntT maxP=50);
    //: Design a shape model given some data.
    //!param: sample    - list of appearances for training the model.
    //!param: variation - percentage of variation preserved during PCA.
    //!param: maxP      - limit on number of parameters contained in the shape model.

    UIntT Dimensions() const
    { return shapeModel.OutputSize() + NoFixedParameters(); }
    //: Size of the parameter vector.

    const SArray1dC<Point2dC> &MeanPoints() const
    { return meanPoints; }
    //: Return mean control point positions (in normalised frame) for shape model.

    SArray1dC<Point2dC> Synthesize(const VectorC &parm) const;
    //: Synthesis a control point set from a parameter vector.

    const VectorC &EigenValues() const
    { return eigenValues; }
    //: Access eigen values.

    const VectorC &FixedMean() const
    { return fixedMean; }
    //: Return mean pose parameters for the shape model.

    void MakePlausible(VectorC &parm, RealT NbSigma = 3) const;
    //: Make 'parm' a plausible parameter vector.
    //  This imposes hard limits of +/-3 std to each parameter.

    bool P2PError(const VectorC &parm,const SArray1dC<Point2dC> &points,VectorC &errVec) const;
    //: Return vector of point to point errors between shape represented by vector 'parm' and target shape defined by 'points'.

	const SampleVectorC GetNormalisedPoints(void) const {return normedpts;}
  protected:
    bool m_verbose; // Write error messages to standard out.

    UIntT nPoints; // Number of control points in the model.

    FunctionC shapeModel;   // Shape model, map control point location to parameters.
    FunctionC invShapeModel;// Inverse shape model, map parameters to control point location

    SArray1dC<Point2dC> meanPoints; // Mean control point positions (in normalised frame).
    VectorC fixedMean; // Mean pose.
	SampleVectorC normedpts;//Projected Points
    VectorC eigenValues;   // eigenValues from PCA.
  };

  //! userlevel=Normal
  //: Statistical shape model.

  class ASMNormalisationC
    : public RCHandleVC<ASMNormalisationBodyC>
  {
  public:
    ASMNormalisationC()
    {}
    //: Default constructor.
    // Creates an invalid handle.

    ASMNormalisationC(BinIStreamC &is);
    //: Binary stream constructor.

    ASMNormalisationC(istream &is);
    //: Stream constructor.

    ASMNormalisationC(bool)
      : RCHandleVC<ASMNormalisationBodyC>(*new ASMNormalisationBodyC())
    {}
    //: Construct a shape model.

  protected:
    ASMNormalisationC(ASMNormalisationBodyC &bod)
      : RCHandleVC<ASMNormalisationBodyC>(bod)
    {}
    //: Body constructor.

    ASMNormalisationC(ASMNormalisationBodyC *bod)
      : RCHandleVC<ASMNormalisationBodyC>(*bod)
    {}
    //: Body ptr constructor.

    ASMNormalisationBodyC &Body()
    { return RCHandleVC<ASMNormalisationBodyC>::Body(); }
    //: Access body.

    const ASMNormalisationBodyC &Body() const
    { return RCHandleVC<ASMNormalisationBodyC>::Body(); }
    //: Access body.

  public:

    void SetVerbose(bool verboseMode)
    { Body().SetVerbose(verboseMode); }
    //: If set to true designer will print error messages to std::cerr

    bool RawParameters(const AAMAppearanceC &inst,VectorC &fixedParams,VectorC &freeParams) const
    { return Body().RawParameters(inst,fixedParams,freeParams); }
    //: Generate raw parameters.
    //!param: inst        - input appearance for which we would like to compute the parameters.
    //!param: fixedParams - output pose parameters (e.g. pose, scale, orientation).
    //!param: freeParams  - output normalised control point coordinates. This vector consists of the concatenation of the X and Y coordinates of all control points in a normalised frame.
    //  The raw parameters are the parameters representing the shape before compression using PCA. They consists of the pose parameters, which describe the pose of the model instance in the image, and the shape parameters, which describe its shape.

    VectorC Parameters(const AAMAppearanceC &inst) const
    { return Body().Parameters(inst); }
    //: Return a parameter vector representing the appearance 'inst'.

    bool Design(const SampleC<AAMAppearanceC> &sample,RealT variation = 0.95,UIntT maxP=25)
    { return Body().Design(sample,variation,maxP); }
    //: Design a shape model given some data.
    //!param: sample    - list of appearances for training the model.
    //!param: variation - percentage of variation preserved during PCA.
    //!param: maxP      - limit on number of parameters contained in the shape model.

    bool Design(SampleStreamC<AAMAppearanceC> &sample,RealT variation = 0.95,UIntT maxP=25)
    { return Body().Design(sample,variation,maxP); }
    //: Design a shape model given some data.
    //!param: sample    - list of appearances for training the model.
    //!param: variation - percentage of variation preserved during PCA.
    //!param: maxP      - limit on number of parameters contained in the shape model.

    UIntT Dimensions() const
    { return Body().Dimensions(); }
    //: Size of the parameter vector.

    SArray1dC<Point2dC> MeanPoints() const
    { return Body().MeanPoints(); }
    //: Return mean control point positions (in normalised frame) for shape model.

    SArray1dC<Point2dC> Synthesize(const VectorC &parm) const
    { return Body().Synthesize(parm); }
    //: Synthesis a control point set from a parameter vector.

    const VectorC &EigenValues() const
    { return Body().EigenValues(); }
    //: Access eigen values.

    IntT NoFixedParameters() const
    { return Body().NoFixedParameters(); }
    //: Return number of parameters describing the pose
    //  These parameters include e.g. the position, scale and orientation of the model instance

    const VectorC &FixedMean() const
    { return Body().FixedMean(); }
    //: Return mean pose parameters for the shape model.

	const SampleVectorC GetNormalisedPoints(void) const
	{ return Body().GetNormalisedPoints();}
	//: Return projected parameters
	
    void MakePlausible(VectorC &parm, RealT NbSigma = 3) const
    { return Body().MakePlausible(parm,NbSigma); }
    //: Make 'parm' a plausible parameter vector.
    //  This imposes hard limits of +/-3 std to each parameter.

    bool P2PError(const VectorC &parm,const SArray1dC<Point2dC> &points,VectorC &errVec) const
    { return Body().P2PError(parm, points, errVec);}
    //: Return vector of point to point errors between shape represented by vector 'parm' and target shape defined by 'points'.
    //: Returns point to point error between shape modelled by vector parm and target shape defined by points
  };

  inline
  BinOStreamC &operator<<(BinOStreamC &s,const ASMNormalisationC &ap) {
    ap.Save(s);
    return s;
  }
  //: Save shape model to binary stream.

  inline
  BinIStreamC &operator>>(BinIStreamC &s,ASMNormalisationC &ap) {
    ap = ASMNormalisationC(s);
    return s;
  }
  //: Read shape model from binary stream.




#endif
