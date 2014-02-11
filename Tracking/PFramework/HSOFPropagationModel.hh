#ifndef HSOPTICALFLOWPROPAGATION_HH
#define HSOPTICALFLOWPROPAGATION_HH
/////////////////////
// File: HSOpticalFlowPropagationModelC
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

class HSOpticalFlowPropagationModelC : public BasePropagationModelC
{
public:
	HSOpticalFlowPropagationModelC(){gaussianwindow = 1;}
	~ HSOpticalFlowPropagationModelC(){}
	HSOpticalFlowPropagationModelC(const PairC<ImageC<RealT> > &ims, const RealT &gw)
	{
		HSOpticFlowC oprop;
		gaussianwindow = gw;
		cout<<"About to get optical flow estimate"<<endl;
		mot_img = oprop.Estimate(ims).Copy();
		cout<<"About to get optical flow estimate"<<endl;
		mean_mot = ComputeMeanMotion(mot_img); //can change this if we need Weighted Mean instead of just mean
		cout<<"Mean Motion = "<<mean_mot<<endl;
	}
	HSOpticalFlowPropagationModelC(const PairC<ImageC<RealT> > &ims)
	{
		HSOpticFlowC oprop;
		gaussianwindow = 1; //default value of the gaussian noise window
		mot_img = oprop.Estimate(ims).Copy();
		mean_mot = ComputeMeanMotion(mot_img); //can change this if we need Weighted Mean instead of just mean
	}	
	// Perform Gaussian Convolution inside class using RealT images
	//Constructor is for application to a pair of filtered images already (so perform Gaussian Convolution
	ParticleC Propagate(ParticleC &pt);
	//Accessor Methods
	ImageC<Vector2dC> GetMotionImage(void) const {return mot_img;}
	Vector2dC GetMeanMotionVector(void) const {return mean_mot;}
	RealT GetNoiseWindow(void) const {return gaussianwindow;}
protected:
	ImageC<Vector2dC> mot_img;
	Vector2dC mean_mot;
	RealT gaussianwindow;
	Vector2dC ComputeMeanMotion(const ImageC<Vector2dC> &im);
	Vector2dC ComputeWeightedMeanMotion(const ImageC<Vector2dC> &im);
};

#endif
