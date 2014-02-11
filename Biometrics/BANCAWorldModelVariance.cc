//This program aims to unbundle the spl biometrics lip feature data created by URS into individual TSV files

//INCLUDE FILES
#include <iostream>
#include <fstream>
#include 	"Ravl/Option.hh"
#include 	"Ravl/OS/Filename.hh"
#include 	"Ravl/OS/Directory.hh"
#include 	"Ravl/PatternRec/Sample.hh"
#include 	"Ravl/PatternRec/SampleIter.hh"
#include 	"Ravl/Vector.hh"
#include 	"Ravl/SumsNd2.hh"
using namespace RavlN;
using namespace std;
//MAIN PROGRAM

SampleC<VectorC> ImportTSVData(const FilenameC& file);

int main(int nargs, char** nargv)
{	
	OptionC opt(nargs,nargv);
	FilenameC worldlist = opt.String("f","BANCAUBM.list","Example world list to input");
	DirectoryC worlddir = opt.String("d","/vol/vssp/lip-tracking/BiometricsWorldModel/","World Model Directory");
	opt.Check();
	
	SampleC<FilenameC> files;
	SumsNd2C mc(22);
	ifstream bancalist;
	bancalist.open(worldlist,ifstream::in);
	if(bancalist.is_open())
	{
		while(! bancalist.eof())
		{
			StringC filename, dummy;
			//filename << bancalist.get();
			bancalist >> filename;
			bancalist >> dummy;
			files.Append((FilenameC)filename.Copy());
		}
	}
	//cout<<files<<endl;
	//Read data in from files
	for(SampleIterC<FilenameC> it(files); it; it++)
	{
		StringC bfile = (StringC)worlddir + (StringC)(*it) + ".TSV";
		SampleC<VectorC> sample = ImportTSVData((FilenameC)bfile);
		for(SampleIterC<VectorC> it2(sample); it2; it2++)
			mc += (*it2);
	}
	cout<<"Mean of Sample is: "<<mc.MeanCovariance().Mean()<<endl;	
	cout<<"Covariance of Sample is: "<<mc.MeanCovariance().Covariance()<<endl;
	return 0;
}
SampleC<VectorC> ImportTSVData(const FilenameC& file)
{
	//Read in file name and import data
	ifstream myfile;
	myfile.open(file.chars(), ifstream::in);
	SampleC<VectorC> res;
	if(myfile.is_open())
	{
		
		UIntT numsamps(0), size(0);
		myfile >> numsamps;
		myfile >> size ;
		VectorC vect(size); vect.Fill(0.0);
		for(UIntT j = 0; j < numsamps; j++)
		{
			for(UIntT i = 0; i < size; i++)
			{
				myfile >> vect[i] ;	
			}
			res.Append(vect.Copy());
		}
		myfile.close();
	}
	return res;
}
