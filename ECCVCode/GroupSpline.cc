#include "GroupSpline.hh"

namespace RavlN
{
void GroupBSpline(Array1dC<ParticleC> &myparticles, SampleC<SArray1dC<Point2dC> > &mycurvepts, SampleC<SArray1dC<LinePP2dC> > &normals)
{
	//Take the input particles and obtain their state vectors from them
	SampleC<VectorC> mycontrolpoints;
	for(Array1dIterC<ParticleC> it(myparticles); it; it++)
	{
		VectorC myvect = PointsToVector(it->GetState().GetObservationVector());
		mycontrolpoints.Append(myvect.Copy());
	}
	//Write the point to a file
	StringC myfile = "/dev/shm/GroupPoints.txt";
	ofstream mystr;
	mystr.open((myfile));
	if(mystr.is_open())
	{
		for(SampleIterC<VectorC> it(mycontrolpoints); it; it++)
		{
			UIntT size = it->Size();
			for(UIntT i = 0; i < size; i++)
			{
				mystr<<(*it)[i]<<"\t";
			}
			mystr<<"\n";
		}
	}
	else
	{
		cerr<<"Could not save the file: "<<myfile<<endl;
		exit(1);
	}
	mystr.close();
	//Call the mathematica Kernal and run the program
	StringC cmd = StringC("/vol/vssp/lip-tracking/ExtSW/Mathematica/Scripts/math  -run \"<< /vol/vssp/lip-tracking/MathematicaNB/groupsplinedata.m\"");
	system(cmd);
	//Now Read in the resulting files
	ifstream mycurvefile("/dev/shm/GroupCurvePoints.txt");
	if(mycurvefile.is_open())
	{
		while(!mycurvefile.eof())
		{
			string line;
			std::getline(mycurvefile,line);
			StringC myline(line);
			StringListC mylinelist(myline, "\t"); 
			VectorC myvect(mylinelist.Size());
			for(UIntT i = 0; i < myvect.Size(); i++)	
				myvect[i] = mylinelist.Nth(i).RealValue();
			if(!myvect[0]==0.)
				mycurvepts.Append(VectorToPoints(myvect.Copy()));
		}
	}
	mycurvefile.close();
	ifstream mynormfile("/dev/shm/GroupNormalPoints.txt");
	SampleC<SArray1dC<Point2dC> > mynormvects;
	if(mynormfile.is_open())
	{
		while(!mynormfile.eof())
		{
			string line;
			std::getline(mynormfile,line);
			StringC myline(line);
			StringListC mylinelist(myline, "\t"); 
			VectorC myvect(mylinelist.Size());
			for(UIntT i = 0; i < myvect.Size(); i++)	
				myvect[i] = mylinelist.Nth(i).RealValue();
			if(!myvect[0]==0.)
				mynormvects.Append(VectorToPoints(myvect.Copy()));
		}
	}
	mynormfile.close();
	//Create the normal stuff
	for(UIntT i = 0; i < mynormvects.Size(); i++)
	{
		SArray1dC<LinePP2dC> arr(mycurvepts.Nth(i).Size());
		for(UIntT j = 0; j < arr.Size(); j++)
		{
			arr[j] = LinePP2dC(Point2dC(mycurvepts.Nth(i)[j]),Vector2dC(mynormvects.Nth(i)[j].Row(),mynormvects.Nth(i)[j].Col()));
		}
		normals.Append(arr.Copy());
	}
	//DONE!
}

}	
