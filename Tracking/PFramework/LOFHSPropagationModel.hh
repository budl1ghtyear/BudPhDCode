#ifndef LOFHSPROPAGATION_HH
#define LOFHSPROPAGATION_HH
/////////////////////
// File: LOFHSPropagationModelC (Local Optical Flow Horn Schunk Optical FLow)
// Author: Bud Goswami
// Date: 28th April 2009
// Purpose: Declare the HS based optical flow based propagation model
/////////////////////

#include "BasePropagationModel.hh"
#include "Particle.hh"
#include "HSOpticFlow.hh"
#include "Ravl/Image/GaussConvolve2d.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Stream.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/Matrix.hh"
#include "Ravl/Pair.hh"
#include "Ravl/RandomGauss.hh"
#include "Ravl/Vector2d.hh"

using namespace RavlN;
using namespace RavlImageN;
using namespace vbapAT;

class LOFHSPropagationModelC : public BasePropagationModelC
{
public:
	LOFHSPropagationModelC(){gaussianwindow = 1;}
	~ LOFHSPropagationModelC(){}
	LOFHSPropagationModelC(const PairC<ImageC<RealT> > &ims, const RealT &gw=1, const RealT &window=2)
	{
		HSOpticFlowC oprop;
		gaussianwindow = gw;
		wind = window;
		cout<<"About to get optical flow estimate"<<endl;
		mot_img = oprop.Estimate(ims).Copy();
		cout<<"About to get optical flow estimate"<<endl;
		//mean_mot = ComputeMeanMotion(mot_img); //can change this if we need Weighted Mean instead of just mean
		//cout<<"Mean Motion = "<<mean_mot<<endl;
	}
	// Perform Gaussian Convolution inside class using RealT images
	//Constructor is for application to a pair of filtered images already (so perform Gaussian Convolution
	ParticleC Propagate(ParticleC &pt);
	//Accessor Methods
	ImageC<Vector2dC> GetMotionImage(void) const {return mot_img;}
	RealT GetNoiseWindow(void) const {return gaussianwindow;}
protected:
	ImageC<Vector2dC> mot_img;
	RealT gaussianwindow;
	IntT wind;
	Vector2dC ComputePatchMeanMotion(const ImageRectangleC &im);
	ImageRectangleC GetPatchRectangle(const Index2dC &ind, const UIntT &window)
	{
		//Given a window, this will return an image rectangle C around the specified index2dC object
		return ImageRectangleC (Max(ind.Row()-window,mot_img.Frame().TRow()),Min(ind.Row()+window,mot_img.Frame().BRow()),
		Max(ind.Col()-window,mot_img.Frame().LCol()),Min(ind.Col()+window,mot_img.Frame().RCol()));
	}
};

#endif
