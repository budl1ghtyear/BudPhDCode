#ifndef LMSOFPROPMODEL_HH
#define LMSOFPROPMODEL_HH

#include "BasePropagationModel.hh"
#include "Ravl/Image/LMSOpticFlow.hh"
#include "Ravl/Pair.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Vector2d.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/Stream.hh"
#include "Ravl/IO.hh"

using namespace RavlN;
using namespace RavlImageN;

class LMSOFPropagationModelC : public BasePropagationModelC
{
	public:
	LMSOFPropagationModelC() {}
	~LMSOFPropagationModelC() {}
	
	LMSOFPropagationModelC(const PairC<ImageC<RealT> > &imgs, const IndexC &patch_window)
	{
		of_obj.SetRegionSize((IntT)patch_window.V());
		ImageC<Vector2dC> m_img1 = of_obj.Estimate(imgs);
		Vector2dC zero_vec(0.0,0.0);
		ImageC<Vector2dC> m_img2(imgs.Data1().Frame());
		m_img2.Fill(zero_vec);
		for(Array2dIter2C<Vector2dC, Vector2dC> it(m_img1, m_img2); it; it++)
		{
			it.Data2() = it.Data1();
		}
		mot_img = m_img2.Copy();
		//Get the drawn optical flow field
		ImageC<ByteYUVValueC> draw_img;of_obj.DrawMotion(imgs[0],draw_img);
		if(!Save("@X:DRAWN MOT_IMG",draw_img)) cerr<<"Could not save drawn image"<<endl;

		window = patch_window;
	}
	//:Constructor with default optical flow estimation parameters
	
	ParticleC Propagate(ParticleC &pt);
	//:Propagation Method
	
	LMSOpticFlowC GetOF(void) const {return of_obj;}	
	//:Return the OpticalFlow method to access the functions inside the OF class
	
	ImageC<Vector2dC> GetOFImage(void) const {return mot_img;}
	//:Return the motion image
	
	protected:
	LMSOpticFlowC  of_obj;
	ImageC<Vector2dC> mot_img;
	IndexC window;
	Vector2dC ComputePatchMeanMotion(const ImageRectangleC &im);
	ImageRectangleC GetPatchRectangle(const Index2dC &ind, const IndexC &window)
	{
		//Given a window, this will return an image rectangle C around the specified index2dC object
		return ImageRectangleC (Max(ind.Row()-window,mot_img.Frame().TRow()),Min(ind.Row()+window,mot_img.Frame().BRow()),
		Max(ind.Col()-window,mot_img.Frame().LCol()),Min(ind.Col()+window,mot_img.Frame().RCol()));
	}	
};

#endif
