#ifndef VBAPAT_HSOPTICFLOW_HEADER
#define VBAPAT_HSOPTICFLOW_HEADER 1

#include "Ravl/Image/Image.hh"
#include "Ravl/Vector2d.hh"
#include "Ravl/Matrix2d.hh"
#include "Ravl/Pair.hh"
#include "Ravl/Image/ByteYUVValue.hh"

namespace vbapAT {

  using namespace RavlN;
  using namespace RavlImageN;
  
  class HSOpticFlowC {
  public:
    HSOpticFlowC (bool verbose=false);
    // Constructor
    
    HSOpticFlowC& SetRegionSize(IntT size)
    { region = size; return *this; }
    // Set size of patch over which motion error is minimised
    
    HSOpticFlowC& SetGradientOrder(IntT order)
    { grad_order = order; return *this; }
    // Set order of spatial gradient estimator
    
    ImageC<Vector2dC> Estimate (const ImageC<Vector2dC> &grad,const ImageC<RealT> &dt);
    // operates on spatial and temporal gradient images

    ImageC<Vector2dC> Estimate (const PairC<ImageC<RealT> > &im);
    // Applies method to pair of (filtered) images using first-order differences

    const ImageC<Vector2dC> &Motion () const
    { return motion; }
    // Return field of motion vectors
    
    static void DrawMotion(const ImageC<RealT> &im,const ImageC<Vector2dC> &motion,ImageC<ByteYUVValueC> &op);
    // Create an image where motion vectors are plotted as U/V components of
    // colour image, with 1st original image as Y component
    
    void DrawMotion(const ImageC<RealT> &im,ImageC<ByteYUVValueC> &op)
    { DrawMotion(im,motion,op); }
    // Create an image where motion vectors are plotted as U/V components of
    // colour image, with 1st original image as Y component
    
  private:
    /* params: */
    RealT alpha;	 // Lagrange multiplier
    IntT iter;		 // number of iterations
    IntT region;         // regression is done over region of this size
    bool verbose;
    IntT grad_order;     // order of spatial gradient (1..3)
    
    /* images: */
    ImageC<Vector2dC> motion; // results of motion estimation
    
  };

  
} 
#endif
