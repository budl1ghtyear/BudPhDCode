//This program aims to unbundle the spl biometrics lip feature data created by URS into individual TSV files

//INCLUDE FILES
#include <iostream>
#include <fstream>
#include 	"Ravl/Option.hh"
#include 	"Ravl/OS/Filename.hh"
#include 	"Ravl/OS/Directory.hh"
#include 	"Ravl/DLIter.hh"
#include 	"Ravl/DList.hh"
#include 	"Ravl/SArray1d.hh"
#include	"Ravl/SArray1dIter.hh"
#include 	"Ravl/PatternRec/Sample.hh"
#include 	"Ravl/PatternRec/SampleIter.hh"
#include 	"Ravl/Vector.hh"
#include 	"BANCAGT.hh"
#include 	"UtilityFunctions.hh"
#include 	"MBspline.hh"
using namespace RavlN;
using namespace std;
//MAIN PROGRAM

VectorC GetVectorFromDIDFile(const FilenameC& file);
void CreateTSVFile(const FilenameC& file, const SampleC<VectorC>& vects);
int main(int nargs, char** nargv)
{	
	OptionC opt(nargs,nargv);
	DirectoryC bancadir = opt.String("d","/vol/vssp/lip-tracking/BANCALIPANNOTATION","Example directory containing BANCA annotations");
	opt.Check();
	
	//Need to hardcode directory structure:
	SArray1dC<DirectoryC> dirs(2);
	dirs[0] = bancadir + "/Controlled/annotation/";
	dirs[1] = bancadir + "/Degraded/annotation/"; 
	for(SArray1dIterC<DirectoryC> it(dirs); it ; it ++)
	{
		DListC<StringC> users = (*it).FiltList("*"); 
		//cout<<"Inside "<<(*it)<<" the program found the following users:\n "<<users<<endl;
		//Now do processing of the internal file:
		//UIntT counter = 0;
		for(DLIterC<StringC> it2(users); it2; it2++)
		{
			//if(counter == 0)
			{
				DListC<StringC> files = DirectoryC((StringC)(*it)+(*it2)).FiltList("*.did");
				//cout<<"Files Inside : "<<DirectoryC((StringC)(*it)+(*it2))<<" are : \n"<<files<<endl;
				FilenameC thisfile((*it2)+".TSV");
				StringC filename = (*it) + (*it2) + "/" + thisfile;
				cout<<"Currently writing: "<<filename<<endl;
				//Now go ahead and create the file:
				SampleC<VectorC> lipvects;
				for(DLIterC<StringC> it3(files); it3 ; it3++)
				{
					StringC filestr = (StringC)(*it) + (*it2) + "/" + (*it3);
					FilenameC didfile = FilenameC(filestr);
					//cout<<"DIDfile - "<<didfile<<endl;
					VectorC vect = GetVectorFromDIDFile(didfile);
					//cout<<vect<<endl;
					lipvects.Append(vect.Copy());
				}
				cout<<"Data to be written is: \n"<<lipvects<<endl;
				//At this point, we have all the information that we need inside the lipvects container. Create the TSV file from this
				CreateTSVFile(filename, lipvects);
			}
			//counter ++;
		}
	}
	return 0;
}

VectorC GetVectorFromDIDFile(const FilenameC& file)
{
	//Load the file and get a dlist of points:
	BANCAGT didfile(file);
	DListC<Point2dC> outer(didfile.GetOuterPoints());
	SArray1dC<Point2dC> outerarr(DListToSArray_B(outer));	//Now we have the input points as a vector, we need to get BSpline control points for them
	MBSplineC spl;
	//cout<<"Input array "<<SArrayToArray_B(outerarr)<<endl;
	Array1dC<Point2dC> ctrlarr = spl.FitQuadraticSpline(SArrayToArray_B(outerarr));
	//cout<<"ctrlarr obtained "<<ctrlarr<<endl;
	return PointsToVector(ArrayToSArray_B(ctrlarr.Copy()));
}

void CreateTSVFile(const FilenameC& file, const SampleC<VectorC>& vects)
{
	ofstream myfile;
	myfile.open (file.chars());
	myfile << vects.Size() << "\t" << vects.First().Size() << "\n";
	for(SampleIterC<VectorC> it(vects); it; it++)
	{
		for(UIntT i = 0; i < (*it).Size(); i++)
		{
			myfile<<(*it)[i]<<"\t";
		}
		myfile<<"\n";
	}
	myfile.close();
}
