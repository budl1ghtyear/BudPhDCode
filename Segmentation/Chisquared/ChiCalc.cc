//Required Libraries
#include "Ravl/DP/DataConv.hh"
#include "Ravl/DP/SequenceIO.hh"
#include "Ravl/DP/FileFormatIO.hh"
#include "Ravl/Math.hh"
#include "Ravl/Option.hh"
//Custom Class Definitions
#include "chisquaredistr.hh"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace RavlN;

int main(int nargs,char **argv) 
{
	//CREATE THE REQUIRED INPUT OPTIONS
	OptionC opt(nargs,argv);
	RealT v = opt.Real("v",1,"Degrees of freedom");//specify the number of subsets for nesting
	RealT y = opt.Real("y",0.5,"Chi-squared Error");//specify the number of subsets for nesting
	
	opt.Check();
	//chisquaredistr ch;
	RealT x = (RealT)invchisquaredistribution(v,y);
	cout<<"D.O.F. = "<<v<<"\t Chi-squared Error = "<<y<<"\t Prob x = "<<x<<endl;
	return 0;
}
