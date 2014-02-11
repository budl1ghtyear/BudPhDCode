#include "Ravl/Assert.hh"
#include "Ravl/DP/SequenceIO.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ByteRGBValue.hh"
#include "Ravl/Image/DeinterlaceStream.hh"

using namespace RavlN;
using namespace RavlImageN;

int main (int argc, char* argv[]) 
{
  RavlAlwaysAssertMsg(argc==3, StringC("Synopsis: \n") + argv[0] + " infile outfile");

  DPIPortC<ImageC<ByteRGBValueC> > in;
  RavlAlwaysAssertMsg(OpenISequence(in, argv[1], "", true),
                      StringC("Could not open ") + argv[1] + " for i/p");

//!Deinterlaced!
  DeinterlaceStreamC<ByteRGBValueC> din(in);

  DPOPortC<ImageC<ByteRGBValueC> > out;
  RavlAlwaysAssertMsg(OpenOSequence(out, argv[2], "", true),
                      StringC("Could not open ") + argv[2] + " for o/p");
  ImageC<ByteRGBValueC> im;
  UIntT i=1;
  cout<<"Start frame: "<<din.Start()<<endl;
  while(din.GetAt(i,im)){
//  while(in.Get(im)){
      out.Put(im);
      cout<<"frame idx: "<<i<<" "<<din.Tell64()<<endl;
      i+=2;
  }
  cout<<"Size: "<<din.Size()/2<<endl;
}
