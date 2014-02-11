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
#include "Ravl/StrStream.hh"
#include "Ravl/Image/DrawFrame.hh"
#include "BoxStateVector.hh"
#include "HistMeasurementModel.hh"
#include "BoxPropModel.hh"
using namespace RavlN;
using namespace RavlImageN;
using namespace std;
RandomGaussC BasePropagationModelC::rnd;

void DrawCrosses(ImageC<RealRGBValueC> &im, SArray1dC<Point2dC> &pts, RealRGBValueC col = RealRGBValueC(255,100,100));
SampleC<Tuple2C<StringC,Index2dC> > LoadBANCAWindow(FilenameC const &file);
Index2dC GetFirstCenter(FilenameC const &v_file, FilenameC const &t_file);

int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	FilenameC vid_file = opt.String("i","./videofile.avi","Input Video File");
	FilenameC test_file = opt.String("f","./firstfile.txt","The coordinates of the mouth rectangle in the first frame");
	//Particle Filter Parameters
	UIntT noise = opt.Int("n",2,"Length of noise perturbtaion in pixels");
	IntT numpts = opt.Int("pt",50,"Number of Particles");
	opt.Compulsory("i");
	opt.Check();
	
	DirectoryC resdir("/vol/vssp/lip-tracking/ThesisTrackExp/BANCABOX/");
	StringC outfilestr = (StringC)resdir+ vid_file.BaseNameComponent()+"_BOX"+StringC(numpts)+StringC(noise)+".txt";
	
	DPIPortC<ImageC<RealRGBValueC> > in;
	cout<<"Trying to open: "<<vid_file<<endl;
	RavlAlwaysAssertMsg(OpenISequence(in, vid_file,"",true),StringC("Could not open ") + vid_file + " for i/p"); 
	cout<<"Loaded video stream"<<endl;
	DeinterlaceStreamC<RealRGBValueC> din(in);
	ImageC<RealRGBValueC> im;
	cout<<"Loaded deinterlaced stream"<<endl;
	UIntT i=1;
	SampleC<Array1dC<Point2dC> > results;
	Array1dC<Point2dC> first_pts; //to store the previous lip co-ordinates
	Array1dC<ParticleC> ptcol;
	while(din.GetAt(i,im))
	{
		//~ if(!Save("@X:Img File", im)) exit(1);
		if(i == 1)
		{
			//If first image, we know the location of the lip points, so we just load it from file
			//Code to locate the centre of first markup
			Index2dC firstcent = GetFirstCenter(vid_file,test_file);
			VectorC st(2);
			st[0] = firstcent.Row();st[1] = firstcent.Col();
			BoxStateVectorC initpts(firstcent);
			RealT initwt = 1.0 / (RealT) numpts;
			ParticleC initparticle(initpts,initwt);
			Array1dC<ParticleC> ptcoln(numpts);
			cout<<"Particle array size: "<<numpts<<endl;
			ptcoln.Fill(initparticle);	
			ptcol = ptcoln.Copy();
			Array1dC<Point2dC> obs = initpts.GetObservationVector();
			//~ DrawFrame(im,RealRGBValueC(255,100,200),IndexRange2dC(obs[0].Row(),obs[3].Row(),obs[0].Col(),obs[3].Col()));    
			//~ if(!Save("@X:MImg File", im)) exit(1);
			results.Append(obs.Copy());
			first_pts = obs.Copy();
		}
		else
		{
			//Do Tracking
			BoxPropModelC ofprop(noise);
			//Crop the image for the search window
			ImageRectangleC myrect = ObtainMouthLocalisationWindow(ArrayToSArray_B(first_pts),1);
			cout<<"INPUT IMAGE RECTANGLE IS: "<<myrect<<endl;
			//ImageC<RealRGBValueC> edgeim(im,myrect);
			//EdgeMeasurementModelC meas(im,myrect,search);
			HistMeasurementModelC meas(im,myrect);
			SystematicResamplingC sysres;
			SIR pfilter(numpts, ofprop, meas, sysres,ptcol);
			first_pts =  pfilter.Apply().GetObservationVector();
			//~ IndexRange2dC todraw(first_pts[0].Row(), first_pts[3].Row(),first_pts[0].Col(), first_pts[3].Col());
			//~ DrawFrame(im,RealRGBValueC(255,100,200),todraw);    
			//~ if(!Save("@X:MImg File", im)) exit(1);
			//DrawCrosses(im, mycurvepoints);
			results.Append(first_pts);
		}
		//~ 
		//~ //if(!Save("@X:First Lips with Img File", im)) exit(1);
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
SampleC<Tuple2C<StringC,Index2dC> > LoadBANCAWindow(FilenameC const &file)
{
	IStreamC mystr(file); 
	SampleC<Tuple2C<StringC,Index2dC> > res;
	while(!mystr.IsEndOfStream())	
	{
		StringC filename;
		mystr >> filename;
		Point2dC tl,br;
		StringC buffer;
		mystr >> tl.Col();
		mystr >> tl.Row();
		mystr >> br.Col();
		mystr >> br.Row();
		Index2dC center((tl.Row() + br.Row()) / 2.,(tl.Col() + br.Col()) / 2.);
		res.Append(Tuple2C<StringC,Index2dC>(filename,center.Copy()));
	}
	mystr.Close();
	return res;
}

Index2dC GetFirstCenter(FilenameC const &v_file, FilenameC const &t_file)
{
	SampleC<Tuple2C<StringC,Index2dC> > markups = LoadBANCAWindow(t_file);
	StringC fname = v_file.BaseNameComponent() + ".ppm";
	Index2dC center;
	for(SampleIterC<Tuple2C<StringC,Index2dC> > it(markups); it; it++)
	{
		if(it->Data1() == fname)
		{
			center = it->Data2().Copy();
			cout<<"Markup found "<<*it<<endl;
			break;
		}
	}
	return center;
}
