//This program aims to unbundle the spl biometrics lip feature data created by URS into individual TSV files

//INCLUDE FILES
#include <iostream>
#include <fstream>
#include 	"Ravl/Option.hh"
#include 	"Ravl/OS/Filename.hh"
#include 	"Ravl/OS/Directory.hh"
using namespace RavlN;
using namespace std;
//MAIN PROGRAM

VectorC GetVectorFromSPL(const FilenameC& file, const UIntT& size=11);

int main(int nargs, char** nargv)
{	
	OptionC opt;
	FilenameC splfile = opt.String("f","000.spl","Example spl file to input")
	UIntT
	opt.Check();
	
	return 0;
}

VectorC GetVectorFromSPL(const FilenameC& file, const UIntT& size)
//This function returns a vector which contains "size" number of 2-D points read in from an spl file
//The vector is going to be in the format:
{
	VectorC vect(size*2);
}
