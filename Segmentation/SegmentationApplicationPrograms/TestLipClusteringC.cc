//LipClusteringC test program

//Just to test the Lip Clustering Test Class
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/Option.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/String.hh"
#include "LipClusteringC.hh"

using namespace RavlN;
using namespace RavlImageN;


int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	FilenameC img_name = opt.String("i","Input.ppm","example input image");
	StringC cluster_type = opt.String("c","CMCD","ClusteringType - values allowed - CMCD,KSMCD,KM,FCM, default = CMCD");
	RealT hvalue = opt.Real("h",0.75,"Confidence Value of MCD - values between 0.5 and 1.0, default = 0.75");
	UIntT numclust = opt.Int("n",3,"Number of clusters for KMeans or Fuzzy C means, default = 3, must be greater than 1");
	opt.Compulsory("i");
	opt.Check();
	
	//Load input image
	ImageC<RealRGBValueC> img;
	if(1==(!Load(img_name, img))) cerr<<"Could not load input image"<<endl;
	//Now we have an input image, lets check what input clustering option was assigned
	UIntT clTypeSel=0;
	Array1dC<StringC> clTypes(4);
	clTypes[0]="CMCD";
	clTypes[1]="KSMCD";
	clTypes[2]="KM";
	clTypes[3]="FCM";
	for (UIntT i=0;i<clTypes.Size();i++)
	{
		if (clTypes[i].matches(cluster_type))
		{
			clTypeSel=i;
			break;
		}
	}
	cout<<"Instantiating class"<<endl;
	LipClusteringC<RealRGBValueC> lc(img,hvalue,numclust);
	cout<<"Running Apply"<<endl;
	DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > trends = lc.Apply((ClusteringType)clTypeSel);	
	cout<<"Found clusters - "<<trends.Size()<<endl;
	return 0;
}
