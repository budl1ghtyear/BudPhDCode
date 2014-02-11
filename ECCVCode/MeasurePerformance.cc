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
#include "Ravl/Polygon2d.hh"
#include "Ravl/Sums1d2.hh"
using namespace RavlN;
using namespace RavlImageN;
using namespace std;
RandomGaussC BasePropagationModelC::rnd;
FilenameC GetDIDFile(FilenameC const &file, UIntT const &num);
RealT GetOverlap(DListC<Point2dC> &gtpts, SArray1dC<Point2dC> &trackpts);
int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	FilenameC res_file = opt.String("i","ResFile.txt","File containing the results of the particle filter experimentation");
	UIntT framediff = opt.Int("n",4,"Number of frames between results");
	opt.Compulsory("i");
	opt.Check();

	//Load the input results file:
	DArray1dC<Tuple2C<SArray1dC<Point2dC>,IntT> > results;
	IStreamC myres(res_file);
	myres >> results;
	SArray1dC<RealT> scores(results.Size()); scores.Fill(0.0);
	Sums1d2C stats;
	//cout<<"Loaded Results: "<<results[1].Data1().Size()<<endl;
	cout<<results.Size()<<endl;
	for(UIntT i = 0; i <= results.Size()-1; i+=framediff)
	{
		//Now we have to get the appropriate filename for the .did file
		FilenameC didf;
		if(framediff == 1)
			didf = GetDIDFile(res_file, i*4);
		else
			didf = GetDIDFile(res_file, i);
		StringListC fname(StringC(res_file.BaseNameComponent()),"_");
		StringC didname = fname.Nth(0)+"_"+fname.Nth(1)+"_"+fname.Nth(2)+"_"+fname.Nth(3)+"_"+fname.Nth(4)+"_"+fname.Nth(5);
		FilenameC didfile = "/vol/vssp/lip-tracking/BANCALIPANNOTATION/Controlled/annotation/"+didname+"/"+StringC(didf);
		//~ cout<<"Didfile: "<<didfile<<endl;
		ImageRectangleC imrec(0.,100.,0.,100.);
		BANCAGT gt(didfile, imrec);
		DListC<Point2dC> pts = gt.GetOuterPoints();
		scores[i] = GetOverlap(pts,results[i].Data1());
		stats+= scores[i];
	}
	cout<<stats.MeanVariance()<<endl;
	return 0;
}
FilenameC GetDIDFile(FilenameC const &file, UIntT const &num)
{
	StringC frame;
	if(num < 10)
		frame = "00"+StringC(num);
	else if(num < 100)
		frame = "0"+StringC(num);
	else
		frame = StringC(num).Copy();
	StringListC fname(StringC(file.BaseNameComponent()),"_");
	StringC didname = fname.Nth(0)+"_"+fname.Nth(1)+"_"+fname.Nth(2)+"_"+fname.Nth(3)+"_"+fname.Nth(4)+"_"+fname.Nth(5)+"_"+frame+".did";
	return FilenameC(didname);
}

RealT GetOverlap(DListC<Point2dC> &gtpts, SArray1dC<Point2dC> &trackpts)
{
	DListC<Point2dC> tpts;
	for(UIntT i =0; i < trackpts.Size(); i++)
	{
		tpts.Append(trackpts[i].Copy());
	}
	Polygon2dC gt(gtpts);
	Polygon2dC est(tpts);
	return gt.Overlap(est);
}
