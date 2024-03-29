#include "ASMRotationNormalisation.hh"
#include "Ravl/PatternRec/SampleIter.hh"
#include "Ravl/VirtualConstructor.hh"
#include "Ravl/DP/FileFormatBinStream.hh"
#include "Ravl/DP/FileFormatStream.hh"
#include "Ravl/TypeName.hh"
#include "Ravl/Moments2d2.hh"
#include "Ravl/Sums1d2.hh"
#include "Ravl/SArray1dIter2.hh"
#include "Ravl/SArray1dIter3.hh"
#include "Ravl/Affine2d.hh"
#include "Ravl/Matrix2d.hh"
#include "Ravl/PatternRec/SampleStreamFromSample.hh"

#define DODEBUG 0
#if DODEBUG
#define ONDEBUG(x) x
#else
#define ONDEBUG(x)
#endif

  void InitScaleRotationShapeModel()
  {}

  //: Default constructor.
  ASMScaleRotationNormalisationBodyC::ASMScaleRotationNormalisationBodyC()
  {}

  //: Load from bin stream.
  ASMScaleRotationNormalisationBodyC::ASMScaleRotationNormalisationBodyC(BinIStreamC &is)
    : ASMNormalisationBodyC(is)
  {}

  //: Load from stream.
  ASMScaleRotationNormalisationBodyC::ASMScaleRotationNormalisationBodyC(istream &is)
    : ASMNormalisationBodyC(is)
  {}

  //: Save to binary stream 'out'.
  bool ASMScaleRotationNormalisationBodyC::Save(BinOStreamC &out) const {
    return ASMNormalisationBodyC::Save(out);
  }

  //: Save to stream 'out'.
  bool ASMScaleRotationNormalisationBodyC::Save(ostream &out) const {
    return ASMNormalisationBodyC::Save(out);
  }

  //: Compute mean control points for the list of appearance provided.
  bool ASMScaleRotationNormalisationBodyC::ComputeMean(SampleStreamC<AAMAppearanceC> &sample) {
    if (sample.Size() == 0)
      return false; // No data points!;

    // Generate initial estimate

    sample.Seek(0);
    AAMAppearanceC appearance;
    sample.Get(appearance);

    nPoints = appearance.Points().Size();
    meanPoints = SArray1dC<Point2dC>(nPoints);

    SArray1dIterC<Point2dC> pit(meanPoints);
    for(pit.First();pit;pit++)
      (*pit) = Point2dC(0,0);

    sample.Seek(0);
    while (sample.Get(appearance)) {
      for(SArray1dIter2C<Point2dC,Point2dC> xit(meanPoints,appearance.Points());xit;xit++)
        xit.Data1() += xit.Data2();
    }
    RealT nSamples = (RealT) sample.Size();
    Moments2d2C moments;
    for(pit.First();pit;pit++) {
      *pit /= nSamples;
      moments += *pit;
    }
    Point2dC mean = moments.Centroid();
    for(pit.First();pit;pit++)
      (*pit) -= mean;


    // Refine means iteratively.

    for(int i = 0;i < 4;i++) {
      SArray1dC<Point2dC> newMeans(nPoints);
      // Set new means to zero.
      for(pit = newMeans;pit;pit++)
        (*pit) = Point2dC(0,0);

      // Got through samples taking mean after rotation correction.

      sample.Seek(0);
      while (sample.Get(appearance)) {
//      for(it.First();it;it++) {
        Affine2dC fit = FitAffine(appearance.Points(),meanPoints);

        Matrix2dC sr = fit.SRMatrix();
        Matrix2dC u,v;
        SVD(sr,u,v);
        Matrix2dC rot = u * v.T(); // Take out scaling.
        Affine2dC norm(rot,fit.Translation());

        for(SArray1dIter2C<Point2dC,Point2dC> xit(newMeans,appearance.Points());xit;xit++)
          xit.Data1() += norm * xit.Data2();
      }

      // Create average
      moments.Reset();
      for(pit.First();pit;pit++) {
        *pit /= nSamples;
        moments += *pit;
      }
      Point2dC mean = moments.Centroid();
      for(pit.First();pit;pit++)
        (*pit) -= mean;

      meanPoints = newMeans;
    }

    return true;
  }

  //: Generate raw parameters.
  //!param: inst        - input appearance for which we would like to compute the parameters.
  //!param: fixedParams - output pose parameters (e.g. pose, scale, orientation).
  //!param: freeParams  - output normalised control point coordinates. This vector consists of the concatenation of the X and Y coordinates of all control points in a normalised frame.
  //  The raw parameters are the parameters representing the shape before applying PCA. They consists of the pose parameters, which describe the pose of the model instance in the image, and the shape parameters, which describe its shape.
  bool ASMScaleRotationNormalisationBodyC::RawParameters(const AAMAppearanceC &inst,VectorC &fixedParams,VectorC &freeParams) const {
    IntT size = inst.Points().Size() * 2;
    freeParams = VectorC (size);
    fixedParams = VectorC(NoFixedParameters());
    SArray1dIterC<Point2dC> pi(inst.Points());

    Moments2d2C moments;
    for(;pi;pi++)
      moments += *pi;

    Point2dC mean = moments.Centroid();
    RealT scale = Sqrt(moments.VarX() + moments.VarY());

    // Sort out parameters with fixed meanings.

    Affine2dC fit = FitAffine(inst.Points(),meanPoints);

    // translation tx, ty
    fixedParams[0] = mean[0];
    fixedParams[1] = mean[1];

    Matrix2dC sr = fit.SRMatrix();

    Matrix2dC u,v;
    SVD(sr,u,v);
    Matrix2dC rot = u * v.T();
    RealT angle = ATan2(rot[0][1],rot[0][0]);

    ONDEBUG(cerr << "Rotation=" << rot << " SR=" << sr << " Det=" << rot.Det() << " Scale=" << scale << " Angle=" << angle << " \n");

    // parameters sx, sy
    // these parameters represent rotation and scaling, the parametrisation
    // is chosen for linearity required by the active search algorithm
    // which does a multivariate linear regression between residual and displacement
    fixedParams[2] = scale * cos(angle) - 1;
    fixedParams[3] = scale * sin(angle);

    rot /= scale;
    SArray1dIterC<RealT> vi(freeParams);

    for(pi.First();pi;pi++) {
      Point2dC p =  (rot * ((*pi) - mean));
      *vi = p[0];
      vi++;
      *vi = p[1];
      vi++;
    }
    return true;

  }

  //: Generate control points defining an appearance from the raw parameters.
  //!param: fixedParams - input pose parameters (e.g. pose, scale, orientation).
  //!param: freeParams  - input normalised control point coordinates. This vector consists of the concatenation of the X and Y coordinates of all control points in a normalised frame.
  //!param: out         - ouput control points
  void ASMScaleRotationNormalisationBodyC::RawProject(const VectorC &fixedParams,const VectorC &freeParams,SArray1dC<Point2dC> &out) const {
    RavlAssert((freeParams.Size()/2) == nPoints);
    if(nPoints != out.Size())
      out = SArray1dC<Point2dC>(nPoints);
    SArray1dIterC<RealT> vi(freeParams);
    Point2dC mean(fixedParams[0],fixedParams[1]);
    Matrix2dC sr(1+fixedParams[2],-fixedParams[3],fixedParams[3],1+fixedParams[2]);

    for(SArray1dIterC<Point2dC> pi(out);pi;pi++) {
      Point2dC p;
      p[0] = (*vi); vi++;
      p[1] = (*vi); vi++;
      (*pi) = (sr * p) + mean;
    }
  }

  //: Return number of parameters describing the pose
  //  These parameters include e.g. the position, scale and orientation of the model instance
  IntT ASMScaleRotationNormalisationBodyC::NoFixedParameters() const {
    return 4;
  }

  RAVL_INITVIRTUALCONSTRUCTOR_FULL(ASMScaleRotationNormalisationBodyC,ASMScaleRotationNormalisationC,ASMNormalisationC);


