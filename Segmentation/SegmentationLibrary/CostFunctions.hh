//CostFunctions.hh
//Author - Bud Goswami
//Purpose - to create a header file that contains the function definitions of region based cost functions using RAVL
#ifndef SEGCostFunctions_HH
#define SEGCostFunctions_HH

//DATA STRUCTURES
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter2.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Image/Image.hh"
#include "SampleData.hh"
#include "Ravl/TFVector.hh"
#include "Ravl/Tuple2.hh"
#include "Ravl/Tuple3.hh"
#include "Ravl/Vector.hh"
#include "Ravl/Vector2d.hh"
//IMAGE STUFF
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/RealHSVValue.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Image/RealYUVValue.hh"
#include "Ravl/IO.hh"
#include "Ravl/Option.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/StdMath.hh"
#include "Ravl/Matrix.hh"
#include "Ravl/Math.hh"
//PATTERN RECOGNITION STUFF
#include "Ravl/MeanCovariance.hh"
#include "Ravl/Image/ConnectedComponents.hh"
#include "Ravl/Image/Segmentation.hh"
//GEOMETRY STUFF
#include "Ravl/Ellipse2d.hh"


using namespace RavlImageN;
using namespace RavlConstN;
//Compute Elliptical Spatial Clusters - Density CG
template <typename T>
DListC<Index2dC> DensityCG(DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > &finalpop,const ImageC<T> &img); 

//Use J based Cluster Labelling and Connected Comp.
DListC<Index2dC> JCC(DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > &clusts,const ImageC<RealRGBValueC> &img); 
DListC<Index2dC> JCC(DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > &clusts,const ImageC<VectorC> &img);
DListC<Index2dC> JCC(DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > &clusts,const ImageC<TFVectorC<RealT,3> > &img); 

bool IsInsideEllipse(const Ellipse2dC &ell, const Index2dC &ind);//For DensityCG function
ImageC<UIntT> GetBinaryLipImage(const DListC<Index2dC> &lip_lst, const ImageRectangleC &imrec);
RealT ComputeSegmentationQuality(const ImageC<UIntT> &gimg, const ImageC<UIntT> &limg);//Scoring Program
RealT ComputeJFunction(const MeanCovarianceC &mc, const VectorC vec);

ImageC<UIntT> LabelImage(const ImageC<RealT> &im);
SegmentationC CentralColumnRegions(const SegmentationC &seg_in);


//Function definition for the templated function HAS to be included over here so that it is available during compilation directly...

template <typename T>
DListC<Index2dC> DensityCG(DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > &finalpop,const ImageC<T> &img)
{
	//cout<<"Finished MCD Algorithm. Clusters found - "<<finalpop.Size()<<endl;
	ImageC<UIntT> m_img(img.Frame(),0);
	//Label a binary image with the different clusters
	UIntT label = 0;
	for(DLIterC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > it(finalpop); it; it++)
	{
		for(DLIterC<Tuple2C<VectorC,Index2dC> > it2((*it).Data2()); it2; it2++)
		{
			m_img[(*it2).Data2()] = label;
		}
		label++;
	}
	SegmentationC dummyseg(m_img, label);
	if(!Save("@X: Test Segmentation Image",dummyseg.RandomImage()))
		cerr<<"Failed to show segmentation image"<<endl;
	
	//Perform Connected Components Analysis on this
	ConnectedComponentsC<UIntT> connected;
	Tuple2C<ImageC<UIntT>, UIntT> result = connected.Apply(m_img); 	
	//Now we have a segmentation image with all zeros except for the lip and noise areas
	SegmentationC seg(result.Data1(),result.Data2());
	//Label Skin as being zero
	for(DLIterC<Tuple2C<VectorC, Index2dC> > it(finalpop.First().Data2()); it; it++)
	{
		seg.SegMap()[(*it).Data2()] = 0;
	}
	seg.RemoveSmallComponents(100);
	if(!Save("@X: Segmentation Image",seg.RandomImage()))
		cerr<<"Failed to show segmentation image"<<endl;
	//Now we have a situation where we can compute the eccentricity of non-zero clusters
	Array1dC<DListC<Index2dC> > clusters(seg.Labels());
	for(Array2dIterC<UIntT> it(seg.SegMap()); it; it++)
	{
		clusters[*it].Append(it.Index().Copy());
	}
	RealT cost = 0.0;
	Index2dC centre = img.Frame().Center();
	DListC<Index2dC> lipind;
	for(Array1dIterC<DListC<Index2dC> > it(clusters); it; it++)
	{
		SArray1dC<Point2dC> elldat((*it).Size());
		IndexC i = 0;
		for(DLIterC<Index2dC> dt((*it)); dt; dt++)
		{
			elldat[i] = (Point2dC)(*dt);
			i++;
		}
		Ellipse2dC ell;
		FitEllipse(elldat,ell);
		//Weight = mass (points inside the ellipse)/ area of the cluster
		UIntT mass = 0;
		for(SArray1dIterC<Point2dC> s(elldat); s; s++)
		{
			Index2dC ellind((*s).Row(), (*s).Col());
			if(IsInsideEllipse(ell, ellind))
				mass++;
		}
		RealT f = 0.0;
		if((mass > 0.0)&&(it.Index() != 0))
			f = (RealT)mass / (RealT)(elldat.Size()*centre.SqrEuclidDistance((Index2dC)ell.Centre())) ; 
		if(f>cost)
		{
			cost = f;
			lipind = (*it).Copy(); 
		}			
	}		
	return lipind;
}
//Perform Bayesian Pixel Labelling Using the 2 largest clusters
//Then perform connected components analysis and return the 2nd largest connected component
template <typename T> 
ImageC<RealT> ClusterGroupingProbSpat(DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > &finalpop,const ImageC<T> &img)
{
	//Extract the chromatic trends information we want
	MeanCovarianceC skin_mc = finalpop.Nth(0);
	MeanCovarianceC lip_mc = finalpop.Nth(1);
	ImageC<RealT> pix_lbl = ProbabilisticPixelLabelling(img,skin_mc,lip_mc);
	return pix_lbl;	
}

//Function to peform Probabilistic Pixel Labelling
template <typename T>
ImageC<RealT> ProbabilisticPixelLabelling(const ImageC<T> &img, const MeanCovarianceC &sk, const MeanCovarianceC &lp)
{
	ImageC<RealT> out_img(img.Frame());
	out_img.Fill(0.0);
	UIntT dim = img[img.Frame().TopLeft()].Size();
	for(Array2dIter2C<T,RealT> it(img, out_img); it; it++)
	{
		//Make a vector object
		VectorC pix_vec(dim);
		for(UIntT i = 0; i < dim; i++)
		{
			pix_vec[i] = it.Data1()[i];
		}
		it.Data2() = ComputeJFunction(lp,pix_vec) - ComputeJFunction(sk,pix_vec);	
	}
	//if(!Save("@X: PPL Image",out_img)) cerr<<"Could not save PPL image"<<endl;
	return out_img;	
}

#endif
