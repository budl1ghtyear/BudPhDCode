//This is the zeroth order PCAStateSpace tracker that uses the active contour measurement model

//INCLUDE FILES:
#include "Ravl/IO.hh"
#include "Ravl/Option.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/DP/SequenceIO.hh"
#include "Ravl/Image/DeinterlaceStream.hh"
#include "Ravl/Assert.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Image/ByteRGBValue.hh"
#include "Ravl/OS/Filename.hh"

using namespace RavlN;
using namespace RavlImageN;

int main(int nargs,char *argv[]) 
{
	OptionC opt(nargs,argv);
	FilenameC vid_file = opt.String("i","./videofile.avi","Input Video File");
	opt.Compulsory("i");
	opt.Check();

	DPIPortC<ImageC<ByteRGBValueC> > in;
	cout<<"Trying to open: "<<vid_file<<endl;
	RavlAlwaysAssertMsg(OpenISequence(in, vid_file,"",true),StringC("Could not open ") + vid_file + " for i/p"); 
	DeinterlaceStreamC<ByteRGBValueC> din(in);
	cout<<"Loaded video stream"<<endl;
	ImageC<ByteRGBValueC> im;
	UIntT i = 1;
	//while(in.Get(im))
	while(din.GetAt(i,im))
	{
		if(!Save("@X:Img File", im)) exit(1);

	}
	return 0;
}
