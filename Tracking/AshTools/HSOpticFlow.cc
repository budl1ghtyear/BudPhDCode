// This file is part of RAVL, Recognition And Vision Library 
// Copyright (C) 2003, University of Surrey
// This code may be redistributed under the terms of the GNU Lesser
// General Public License (LGPL). See the lgpl.licence file for details or
// see http://www.gnu.org/copyleft/lesser.html
// file = HSOpticFlow.cc
// The source is based on Horn & Schunck's method.
// Created on 11/04/2009 by Ashish Doshi.
// Modified on 27/04/2009.
// file-header-ends-here

// global include
#include "Ravl/StdConst.hh"
#include "Ravl/Vector2d.hh"
#include "Ravl/Matrix2d.hh"
#include "Ravl/VectorMatrix2d.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/Array2dIter4.hh"
#include "Ravl/Array2dIter2.hh"
#include "Ravl/Array2dIter3.hh"
#include "Ravl/Array2dIter5.hh"
#include "Ravl/Array2dIter7.hh"
#include "Ravl/StdMath.hh"
#include "Ravl/Stream.hh"
#include "Ravl/StrStream.hh"
#include "Ravl/DP/FileFormatIO.hh"
#include "Ravl/Pair.hh"
#include "Ravl/Image/SumRectangles.hh"
#include "Ravl/Image/ScaleValues.hh"
#include "Ravl/Image/SpatialDifference.hh"
#include "Ravl/Image/Convolve2d.hh"

// local include
#include "./HSOpticFlow.hh"

namespace vbapAT {


  // default constructor  
  HSOpticFlowC::HSOpticFlowC (bool Verbose)
    : alpha (0.12),
      iter (5),//changed from 50
      region (3),
      verbose (Verbose),
      grad_order(2)
  {}

  ImageC<Vector2dC> HSOpticFlowC::Estimate (const ImageC<Vector2dC> &grad, const ImageC<RealT> &dt)
  {
    // compute intersection of gradient fields
    ImageRectangleC rect = grad.Frame();
    rect.ClipBy(dt.Frame());
    
    // compute outer products of grad I, dI/dt
    // =======================================
    ImageC<Matrix2dC> grad_grad(rect);
    ImageC<Vector2dC> grad_t (rect);
    ImageC<RealT> dt_sq (rect);
    
    // grad_t = grad * -dt
    //cerr << "    CP 1" << endl;
    for(Array2dIter5C<Vector2dC,Vector2dC,RealT,RealT,Matrix2dC> it(grad_t,    // 1
								    grad,      // 2
								    dt,        // 3
								    dt_sq,     // 4
								    grad_grad, // 5
								    rect);it;it++) {
      it.Data2().OuterProduct(it.Data5());
      Mul(it.Data2(),-it.Data3(),it.Data1());
      it.Data4() = Sqr(it.Data3());
    }
    
    // find the sums of products
    // ======================
    
    //IndexRange2dC mask(Index2dC(0,0),region,region); // Create mask of region by region centered on 0,0 
    
    //ImageC<Matrix2dC> sum_grad_grad;
    //SumRectangles(grad_grad,mask,sum_grad_grad);
    
    //ImageC<Vector2dC> sum_grad_t;
    //SumRectangles(grad_t ,mask,sum_grad_t);
    
    //ImageC<RealT> sum_dt_sq;
    //SumRectangles(dt_sq,mask,sum_dt_sq);
    
    // hence solve for motion
    // =======================================HSOpticFlow.cc=====
    
    //ImageRectangleC sum_rect (sum_grad_grad.Rectangle());
    ImageRectangleC sum_rect (grad.Rectangle());
    if (!sum_rect.IsValid()) {
      if (verbose) cerr << "Images too small or filters too large to compute motion\n";
      motion = ImageC<Vector2dC>();
      return motion;
    }

    //cerr << "    CP 2" << endl;
    motion = ImageC<Vector2dC>(sum_rect);
    ImageC<Vector2dC> mave1,mave2;

    // Create mask for average velocities
    // ============================================
    Array2dC<RealT> k1;
    Array2dC<RealT> k2;
    const RealT mult1=1.0/6.0;
    const RealT mult2=1.0/12.0;
    StrIStreamC ("-1 1 -1 1 0 1 0 1 0 1 0 1 0") >> k1;
    StrIStreamC ("-1 1 -1 1 1 0 1 0 0 0 1 0 1") >> k2;
    //cerr << "    CP 3" << endl;
    Convolve2dC<Vector2dC,Vector2dC,RealT,Vector2dC> conv1(k1);
    Convolve2dC<Vector2dC,Vector2dC,RealT,Vector2dC> conv2(k2);
    
    
    // Iterative solution
    // ============================================
    Vector2dC conv_mean(0,0);Vector2dC sum_v(0.0, 0.0);
    Vector2dC min_diff(0.07,0.07);const IntT s_area = sum_rect.Area();
    for (IntT itt=0; itt<iter; itt++)
    {
		do
		{
        // apply Laplacian filter
        // ========================================
        conv1.Apply(motion,mave1);
	conv2.Apply(motion,mave2);
        //cerr << "    CP 4" << endl;
	mave1 = mave1*Vector2dC(mult1,mult1) + mave2*Vector2dC(mult2,mult2);

        // minimisation
        // ========================================
	
    	RealT sum_v_sq(0.0);
        //cerr << "    CP 5" << endl;
	//cerr << grad.Frame() << endl;
	//cerr << grad_t.Frame() << endl;
	//cerr << mave1.Frame() << endl;
	//cerr << motion.Frame() << endl;
	//cerr << rect.Size() << endl;
        for (Array2dIter4C<Vector2dC,Vector2dC,Vector2dC,Vector2dC> mit(mave1,  	// 1
				 					grad_t,         // 2
                                 					grad,		// 3
				 					motion,		// 4
				 					false);mit;mit++) {
		// compute motion vector
		Vector2dC v;
		const RealT sqr_alpha = Sqr(alpha);
		//cerr << "     CP 5a" << endl;
		RealT denom1 = mit.Data3()[0] / (Sqr(mit.Data3()[0]) + Sqr(mit.Data3()[1]) + sqr_alpha);
		RealT denom2 = mit.Data3()[1] / (Sqr(mit.Data3()[0]) + Sqr(mit.Data3()[1]) + sqr_alpha);
                //cerr << "     CP 5b" << endl;
		v[0] = mit.Data1()[0] - (denom1*(mit.Data3()[0]*mit.Data1()[0] + mit.Data3()[1]*mit.Data1()[1] + mit.Data2()[0]));
		v[1] = mit.Data1()[1] - (denom2*(mit.Data3()[0]*mit.Data1()[0] + mit.Data3()[1]*mit.Data1()[1] + mit.Data2()[1]));
		//cerr << "     CP 5c" << endl;
		mit.Data4() = v;
		sum_v += v;
		sum_v_sq += v.Dot(v);
	}
	if (verbose) {	      
	      cerr << "Mean V = " << sum_v/s_area << '\n';
	      cerr << "S.D. = " << sqrt(sum_v_sq/s_area - Vector2dC(sum_v/s_area).SqrNorm())<< '\n';
	}
	}while((Abs((sum_v/s_area - conv_mean)[0]) > min_diff[0])&&(Abs((sum_v/s_area - conv_mean)[1]) > min_diff[1]));
	cout<<"itt = "<<itt<<"\t Conv_mean = "<<conv_mean<<endl;
	conv_mean = (sum_v/s_area).Copy();
   }   
    //cerr << "    CP 6" << endl;
    return motion;
  }
  
  
  ImageC<Vector2dC> HSOpticFlowC::Estimate (const PairC<ImageC<RealT> > & im)
  {
    // compute spatial image gradient from mean image
    if (((IntT)im[0].Rows() < 2*grad_order+1) || ((IntT)im[0].Cols() < 2*grad_order+1)) {
      cerr << "Image too small for gradient calculation\n";
      return ImageC<Vector2dC>();
    }
    ImageC<RealT> aver(im[0].Frame());
    ImageC<RealT> diff(im[0].Frame());
    for(Array2dIter4C<RealT,RealT,RealT,RealT> it(aver,diff,im[0],im[1]);it;it++) {
      it.Data1() = (it.Data3() + it.Data4())/2;
      it.Data2() = it.Data3() - it.Data4();
    }
    if(!Save("@X:Inside HSOF Aver",aver)) cerr << "Could not show average image" << endl;
	if(!Save("@X:Inside HSOF Diff",diff)) cerr << "Could not show difference image" << endl;
    ImageC<Vector2dC> grad;
    SpatialDifference(grad_order,aver,grad);
    
    // compute motion from spatial gradient and frame difference
    return Estimate(grad, diff);
  }
/*
  
  // Return field of eigenvalues of spatial gradient outer product
  
  ImageC<Vector2dC> LMSOpticFlowC::Eigenvalues() const { 
    ImageC<Vector2dC> ret(lambda.Frame());
    RealT div = (RealT) (region*region);
    for(Array2dIter2C<Vector2dC,Vector2dC> it(ret,lambda);it;it++)
      it.Data1() = it.Data2() / div;
    return ret;
  }
*/  
  void HSOpticFlowC::DrawMotion(const ImageC<RealT> &im,const ImageC<Vector2dC> &motion,ImageC<ByteYUVValueC> &op) {
    // find max motion value
    RealT max(0.0);
    for (Array2dIterC<Vector2dC> i(motion); i.IsElm(); i.Next())
      if (i.Data().SqrNorm()>max)  max = i.Data().SqrNorm();
    RealT factor(127.0/sqrt(max));
    // write motion as colour to o/p image, along with original image in faint
    // grey
    
    ImageC<RealT> scaled;
    ScaleOffsetValues(ImageC<RealT>(im,motion.Frame()),0.0,128.0,scaled);
    
    op = ImageC<ByteYUVValueC>(motion.Rectangle());
    for (Array2dIter3C<ByteYUVValueC,Vector2dC,RealT> it(op, motion, scaled);it;it++)
      it.Data1() = ByteYUVValueC((IntT) it.Data3() + 128,(IntT) (it.Data2()[0]*factor), (IntT) (it.Data2()[1]*factor));
    // add a colour circle key
    
    IndexC keysize(op.Cols()/10);
    if (keysize > op.Rows())  
      keysize = op.Rows();
    IndexRange2dC rect(op.TRow(), op.TRow()+keysize-1, op.RCol()-keysize+1, op.RCol());
    ImageC<ByteYUVValueC> key(op, rect);
    Index2dC centre(rect.Center());
    for (Array2dIterC<ByteYUVValueC> itk(key); itk.IsElm(); itk.Next()) {
      itk.Data().Y() = 255;
      itk.Data().U() = (SByteT) ((itk.Index()-centre)[0]*255.0/keysize);
      itk.Data().V() = (SByteT) ((itk.Index()-centre)[1]*255.0/keysize);
    }
  }

}
