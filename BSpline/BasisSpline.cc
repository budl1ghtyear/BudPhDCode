#include "BasisSpline.hh"

const StringC BasisSplineC::Mathematica = "/vol/vssp/lip-tracking/ExtSW/Mathematica/Scripts/math";
const StringC BasisSplineC::ScriptPath = "/vol/vssp/lip-tracking/MathematicaNB/";
const StringC BasisSplineC::ResPath = "/dev/shm/";
	

BasisSplineC::BasisSplineC(UIntT const &ncp, UIntT const &deg, UIntT const &num):numcp(ncp),degree(deg),numpts(num)
{
	//Create a zero valued control point array
	SArray1dC<Point2dC> pts(ncp);
	Point2dC zero(0.0,0.0); pts.Fill(zero);
	controlpoints = pts.Copy();
}

//Constructor to use when perform curve generation
BasisSplineC::BasisSplineC(SArray1dC<Point2dC> const &ctrlpts, UIntT const &deg, UIntT const &num):controlpoints(ctrlpts.Copy()),numcp(ctrlpts.Size()),degree(deg),numpts(num)
{
	this->GenerateCurve(this->controlpoints);
}

bool BasisSplineC::LoadPointsFromFile(FilenameC const &f_name, UIntT const &size, SArray1dC<Point2dC> &res)
{
	IStreamC myfile(f_name);
	if(myfile.IsOpen())
	{
		SArray1dC<Point2dC> pts(size);
		UIntT index = 0;
		while((!myfile.eof())&&(index < size))
		{
			Point2dC pt(0.0,0.0);
			myfile>>pt.Row();
			myfile>>pt.Col();
			pts[index++] = pt.Copy();
		}
		myfile.Close();
		res = pts.Copy();
		return 1;
	}
	else
	{
		cerr<<"Could not load from the file "<<f_name<<endl;
		return 0;
	}
}

bool BasisSplineC::WritePointsToFile(SArray1dC<Point2dC> const &pts, FilenameC const &f_name)
{
	ofstream myfile;
	myfile.open (f_name);
	if(myfile.is_open())
	{
		for(SArray1dIterC<Point2dC> it(pts); it; it++)
		{
			myfile<<(*it).Row()<<"\t"<<(*it).Col()<<endl;
		}
		myfile.close();
		return 1;
	}	
	else
	{
		cout<<"Could not open the file: "<<f_name<<" for writing"<<endl;
		return 0;
	}
}

bool BasisSplineC::WriteIntToFile(UIntT const &n, FilenameC const &f_name)
{
	ofstream myfile;
	myfile.open(f_name);
	if(!myfile.is_open())
	{
		cerr<<"Could not open the destination file for writing : "<<f_name<<endl;
		return 0;		
	}
	else
	{
		myfile<<n;
		myfile.close();
		return 1;
	}
}

void BasisSplineC::CallSystem(StringC const &cmd)
{
	UIntT i = system(cmd);
	if(i != 0)
	{
		cout<<"Could not execute command - "<<cmd<<endl;
		exit(1);
	}
}

void BasisSplineC::RemoveFiles(FilenameC const &file)
{
	CallSystem("rm "+file);
}

void BasisSplineC::FitCurve(SArray1dC<Point2dC> const &pts)
{
	//Create the necessary writing files:
	FilenameC lipptsfile = ResPath + "000LipPoints.txt";
	FilenameC numcpfile = ResPath + "000Numcp.txt";
	FilenameC numptsfile = ResPath + "000NumPoints.txt";
	FilenameC degreefile = ResPath + "000Degree.txt";
	if(WritePointsToFile(pts,lipptsfile) && WriteIntToFile(this->numcp,numcpfile) && WriteIntToFile(this->numpts,numptsfile) && WriteIntToFile(this->degree,degreefile))
	{
		//run the mathematica script
		StringC cmd = Mathematica + " -run \"<<" + ScriptPath + "SplineFit.m\"";
		//cout<<"About to run the Mathematica Command: "<<cmd<<endl;
		CallSystem(cmd);
		//Create the files to read from:
		FilenameC ctrlptsfile = ResPath + "000ControlPoints.txt";
		FilenameC curveptsfile = ResPath + "000CurvePoints.txt";
		FilenameC normalptsfile = ResPath + "000NormalPoints.txt";
		SArray1dC<Point2dC> norms;
		if(!(LoadPointsFromFile(ctrlptsfile,numcp,controlpoints) && LoadPointsFromFile(curveptsfile,numpts,curvepoints) && LoadPointsFromFile(normalptsfile,numpts,norms)))
		{
			cerr<<"Could not load the files to read in data for B-Spline Curve Fitting"<<endl;
			exit(1);
		}
		else
		{
			//control points and curve points have been instantiated, so now convert the norms into the appropriate structure:
			SArray1dC<LinePP2dC> lines(norms.Size());
			for(SArray1dIter3C<Point2dC, Point2dC, LinePP2dC> it(curvepoints,norms,lines); it; it++)
			{
				Vector2dC vec(it.Data2().Row(),it.Data2().Col());
				LinePP2dC l(it.Data1(),(Vector2dC)vec.Copy());
				it.Data3() = l;
			}
			curvenormals = lines.Copy();
		}
		//now that the entire structure has been instantiated, perform the deletion of the footprint:
		//~ RemoveFiles(lipptsfile);RemoveFiles(numcpfile);RemoveFiles(numptsfile);RemoveFiles(degreefile);
		//~ RemoveFiles(ctrlptsfile);RemoveFiles(curveptsfile);RemoveFiles(normalptsfile);
	}
	else
	{
		cerr<<"Could not write the B-Spline data to text files on the system"<<endl;
		exit(1);
	}
}


void BasisSplineC::GenerateCurve(SArray1dC<Point2dC> const &ctrlpts)
{
		//Create the necessary writing files:
	FilenameC ctrlptsfile = ResPath + "000ControlPoints.txt";
	FilenameC numcpfile = ResPath + "000Numcp.txt";
	FilenameC numptsfile = ResPath + "000NumPoints.txt";
	FilenameC degreefile = ResPath + "000Degree.txt";
	if(WritePointsToFile(ctrlpts,ctrlptsfile) && WriteIntToFile(this->numcp,numcpfile) && WriteIntToFile(this->numpts,numptsfile) && WriteIntToFile(this->degree,degreefile))
	{
		//run the mathematica script
		StringC cmd = Mathematica + " -run \"<<" + ScriptPath + "SplineGenerate.m\"";
		//cout<<"About to run the Mathematica Command: "<<cmd<<endl;
		CallSystem(cmd);
		//Create the files to read from:
		controlpoints = ctrlpts.Copy();
		FilenameC curveptsfile = ResPath + "000CurvePoints.txt";
		FilenameC normalptsfile = ResPath + "000NormalPoints.txt";
		SArray1dC<Point2dC> norms;
		if(!(LoadPointsFromFile(curveptsfile,numpts,curvepoints) && LoadPointsFromFile(normalptsfile,numpts,norms)))
		{
			cerr<<"Could not load the files to read in data for B-Spline Curve Fitting"<<endl;
			exit(1);
		}
		else
		{
			//control points and curve points have been instantiated, so now convert the norms into the appropriate structure:
			SArray1dC<LinePP2dC> lines(norms.Size());
			//cout<<"Normal Line Obtained is : "<<norms<<endl;
			for(SArray1dIter3C<Point2dC, Point2dC, LinePP2dC> it(curvepoints,norms,lines); it; it++)
			{
				Vector2dC vec(it.Data2().Row(),it.Data2().Col());
				LinePP2dC l(it.Data1(),(Vector2dC)vec.Copy());
				it.Data3() = l;
			}
			curvenormals = lines.Copy();
		}
		//now that the entire structure has been instantiated, perform the deletion of the footprint:
		//~ RemoveFiles(numcpfile);RemoveFiles(numptsfile);RemoveFiles(degreefile);
		//~ RemoveFiles(ctrlptsfile);RemoveFiles(curveptsfile);RemoveFiles(normalptsfile);
	}
	else
	{
		cerr<<"Could not write the B-Spline data to text files on the system"<<endl;
		exit(1);
	}
}
