//This is the zeroth order PCAStateSpace tracker that uses the active contour measurement model

//INCLUDE FILES:
#include "Ravl/Option.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/DP/SequenceIO.hh"
#include "Ravl/Image/DeinterlaceStream.hh"
#include "Ravl/PatternRec/Sample.hh"
#include "Ravl/Tuple2.hh"
#include "BasisSpline.hh"
#include "Particle.hh"
#include "SystematicResampling.hh"
#include "SIR.hh"
#include "TrackingUtilityFunctions.hh"
#include "PCAStateVector.hh"
#include "ZOAVPropModel.hh"
#include "ACMeasurementModel.hh"
#include "EdgeMeasurementModel.hh"
#include "ICAStateVector.hh"
#include "BANCAGT.hh"
#include "Ravl/Assert.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Image/DeinterlaceStream.hh"
#include "Ravl/PatternRec/Sample.hh"
#include "Ravl/PatternRec/SampleIter.hh"
#include "Ravl/StringList.hh"
#include <string>
#include "GroupSpline.hh"
using namespace RavlN;
using namespace RavlImageN;
using namespace std;
RandomGaussC BasePropagationModelC::rnd;

void DrawCrosses(ImageC<RealRGBValueC> &im, SArray1dC<Point2dC> &pts, RealRGBValueC col = RealRGBValueC(255,100,100));


int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	FilenameC vid_file = opt.String("i","./videofile.avi","Input Video File");
	FilenameC first_file = opt.String("f","./firstfile.did","The coordinates of the first lip shape that was bootstrapped. This is in .did format");
	DirectoryC pca_dir = opt.String("p","/vol/vssp/lip-tracking/StatisticalShapeModelsOfLips/JADE/RedJADE/","PCA Data Directory");
	UIntT subspace_proj = opt.Int("sp",0,"Method of subspace projection to use- 0: PCA, 1:ICA");
	UIntT numcomponents = opt.Int("npca",10,"Number of eigen-components to choose in state representation");
	//Particle Filter Parameters
	UIntT search = opt.Int("sw",20,"Length of search window in pixels in affine space");
	UIntT noise = opt.Int("n",2,"Length of noise perturbtaion in pixels");
	IntT numpts = opt.Int("pt",50,"Number of Particles");
	opt.Compulsory("i");
	opt.Check();
	
	DirectoryC resdir("/vol/vssp/lip-tracking/ThesisTrackExp/");
	StringC outfilestr = (StringC)resdir+first_file.BaseNameComponent()+"_ICA"+StringC(numpts)+StringC(search)+".txt";
	
	DPIPortC<ImageC<RealRGBValueC> > in;
	cout<<"Trying to open: "<<vid_file<<endl;
	RavlAlwaysAssertMsg(OpenISequence(in, vid_file,"",true),StringC("Could not open ") + vid_file + " for i/p"); 
	cout<<"Loaded video stream"<<endl;
	DeinterlaceStreamC<RealRGBValueC> din(in);
	ImageC<RealRGBValueC> im;
	cout<<"Loaded deinterlaced stream"<<endl;
	UIntT i=1;
	SampleC<Tuple2C<SArray1dC<Point2dC>,IntT> > results;
	SArray1dC<Point2dC> first_pts; //to store the previous lip co-ordinates
	Array1dC<ParticleC> ptcol;
	while(din.GetAt(i,im))
	{
		//if(!Save("@X:Img File", im)) exit(1);
		if(i == 1)
		{
			//If first image, we know the location of the lip points, so we just load it from file
			BANCAGT bancabootstrap(first_file,im.Frame());
			DListC<Point2dC> banca_outerpts = bancabootstrap.GetOuterPoints();
			SArray1dC<Point2dC> dummy(banca_outerpts.Size());
			for(UIntT pt_ind = 0; pt_ind < banca_outerpts.Size(); pt_ind++)
			{
				dummy[pt_ind] = banca_outerpts.Nth(pt_ind).Copy();
			}
			first_pts = dummy.Copy();
			DrawCrosses(im, first_pts);
			cout<<"First_pts: "<<first_pts<<endl;
			//Initialise the system
			BasisSplineC spl(11,2,30);
			spl.FitCurve(first_pts);
			//cout<<"B-Spline Control Points: "<<spl.GetControlPoints()<<endl;
			ICAStateVectorC initpts(pca_dir, spl.GetControlPoints());
			RealT initwt = 1.0 / (RealT) numpts;
			ParticleC initparticle(initpts,initwt);
			Array1dC<ParticleC> ptcoln(numpts);
			cout<<"Particle array size: "<<numpts<<endl;
			ptcoln.Fill(initparticle);	
			ptcol = ptcoln.Copy();
			//~ //StateVector Experiment:
			//cout<<"From inside PCA Vector: "<<initpts.GetObservationVector()<<endl;
			//~ DrawCrosses(im, first_pts);
			//~ SArray1dC<Point2dC> splcp = ArrayToSArray_B(initpts.GetObservationVector());
			//~ BasisSplineC myspl(splcp);
			//~ SArray1dC<Point2dC> mycurvepoints = myspl.GetCurvePoints();
			//~ DrawCrosses(im, mycurvepoints);
			//~ SimilarityProjectionC sim("/vol/vssp/lip-tracking/StatisticalShapeModelsOfLips/PCA/");
			//~ Tuple2C<SArray1dC<Point2dC>,VectorC> affpts = sim.ProjectToAffineSpace(first_pts);
			//~ cout<<"Projection into affine space: "<<affpts<<endl;
			//~ SArray1dC<Point2dC> simres = sim.ProjectFromAffineSpace(affpts.Data1(), affpts.Data2());
			//~ cout<<"Projection from affine space: "<<simres<<endl;
			//~ SampleC<SArray1dC<Point2dC> > curvepts;
			//~ SampleC<SArray1dC<LinePP2dC> > normallines;
			//~ GroupBSpline(ptcoln,curvepts,normallines);
			//~ cout<<"Num Curve Pts "<<curvepts.Size()<<"\t Num normal lines: "<<normallines.Size()<<"\n Normal Lines: \n"<<normallines.Nth(49)<<endl;
			//exit(1);
			results.Append(Tuple2C<SArray1dC<Point2dC>,IntT>(first_pts,i));
		}
		else
		{
			//Do Tracking
			ZOAVPropModelC ofprop(noise);
			//Crop the image for the search window
			ImageRectangleC myrect = ObtainMouthLocalisationWindow(first_pts,search);
			cout<<"INPUT IMAGE RECTANGLE IS: "<<myrect<<endl;
			ImageC<RealRGBValueC> edgeim(im,myrect);
			EdgeMeasurementModelC meas(im,myrect,search);
			SystematicResamplingC sysres;
			SIR pfilter(numpts, ofprop, meas, sysres,ptcol);
			Array1dC<Point2dC> pf_pts = pfilter.ApplyGroup().GetObservationVector(); //GET THE APPLY METHOD 
			first_pts = ArrayToSArray_B(pf_pts);
			BasisSplineC myspl(first_pts);
			SArray1dC<Point2dC> mycurvepoints = myspl.GetCurvePoints();
			//DrawCrosses(im, mycurvepoints);
			results.Append(Tuple2C<SArray1dC<Point2dC>,IntT>(first_pts,i));
		}
		
		//if(!Save("@X:First Lips with Img File", im)) exit(1);
		i+=2;
	}
	OStreamC myres(outfilestr);
	//results>>myres;
	results.DArray().Save(myres);
	return 0;
}

void DrawCrosses(ImageC<RealRGBValueC> &im, SArray1dC<Point2dC> &pts, RealRGBValueC col)
{
	//RealRGBValueC col(255,10,10);
	for(SArray1dIterC<Point2dC> it(pts); it; it++)
	{
		Index2dC ind((*it).Row(),(*it).Col());
		DrawCross(im,col,ind,4);
	}
	//if(!Save("@X:Cross Image",im)) cerr<<"Could not display crosses image"<<endl;
}
