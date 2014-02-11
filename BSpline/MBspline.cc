#include "MBspline.hh"

const StringC MBSplineC::Mathematica = "/vol/vssp/lip-tracking/ExtSW/Mathematica/Scripts/math";
const StringC MBSplineC::ScriptPath = "/vol/vssp/lip-tracking/MathematicaNB/";
const StringC MBSplineC::ResPath = "/dev/shm/LipExp/";
Array1dC<Point2dC> MBSplineC::FitQuadraticSpline(Array1dC<Point2dC> const &pts)
{
	//Write the points to an output file just as space separated points x1 y1 \n x2 y2
	//cout<<"INSIDE FIT QUADRATIC SPLINE METHOD: "<<endl;
	//cout<<pts<<endl;
	FilenameC ptfile(ResPath + "000LipPoints.txt");
	cout<<ptfile<<endl;
	this->WritePointsToFile(pts,ptfile);
	//Create a command file
	StringC cmd = Mathematica + " -run \"<<" + ScriptPath + "make_control_points.m\"";
	CallSystem(cmd);
	//Now load in the input file
	FilenameC resfile = ResPath + "000LipPoints.txt.csv";
	Array1dC<Point2dC> ctrlPts = LoadPointsFromFile(resfile,11);
	//cout<<ctrlPts<<endl;
	//Clean up residuals
	StringC rmcmd1 = "rm "+ptfile; 
	StringC rmcmd2 = "rm "+resfile; 
	CallSystem(rmcmd2);	
	CallSystem(rmcmd1);
	//this->SetControlPoints(ctrlPts);
	//this->SetDegree(2);
	return ctrlPts;
}

Array1dC<Point2dC> MBSplineC::RenderClosedQuadraticSpline(Array1dC<Point2dC> const &pts, UIntT const &npts)
{
	//Write the data points into a file
	FilenameC ptfile = ResPath + "000ContourPoints.txt"; WritePointsToFile(pts,ptfile);
	FilenameC numpt = ResPath + "000NumPoints.txt";
	ofstream myfile; myfile.open(numpt);
	if(myfile.is_open())
		myfile<<npts;
	myfile.close();
	//Now run Mathematica Command
	StringC cmd = Mathematica + " -run \"<<" + ScriptPath + "rendersplinecurve.m\"";
	CallSystem(cmd);
	//Now load in the resulting input file
	FilenameC resfile = ResPath + "000CurvePoints.txt";
	Array1dC<Point2dC> curvePts = LoadPointsFromFile(resfile,npts);
	//Clean up residuals
	CallSystem("rm "+resfile);CallSystem("rm "+ptfile);CallSystem("rm " + numpt);
	return curvePts;
}

Array1dC<LinePP2dC> MBSplineC::GetNormalLines(Array1dC<Point2dC> const &pts, UIntT const &npts)
{
	//Write the data points into a file
	FilenameC ptfile = ResPath + "000ContourPoints.txt"; WritePointsToFile(pts,ptfile);
	FilenameC numpt = ResPath + "000NumPoints.txt";
	ofstream myfile; myfile.open(numpt);
	if(myfile.is_open())
		myfile<<npts;
	myfile.close();
	//Now run Mathematica Command
	StringC cmd = Mathematica + " -run \"<<" + ScriptPath + "getderivative.m\"";
	CallSystem(cmd);
	//Now load in the resulting input files:
	FilenameC resfile = ResPath + "000CurvePoints.txt";
	Array1dC<Point2dC> curvePts = LoadPointsFromFile(resfile,npts);
	FilenameC normfile = ResPath + "000NormalPoints.txt";
	Array1dC<Point2dC> normPts = LoadPointsFromFile(normfile,npts);
	Array1dC<LinePP2dC> lines(curvePts.Size());
	for(Array1dIter3C<Point2dC, Point2dC, LinePP2dC> it(curvePts,normPts,lines); it; it++)
	{
		Vector2dC vec(it.Data2().Row(),it.Data2().Col());
		LinePP2dC l(it.Data1(),(Vector2dC)vec.Copy());
		it.Data3() = l;
	}
	//Clean up residuals
	CallSystem("rm "+resfile);CallSystem("rm "+ptfile);CallSystem("rm " + numpt); CallSystem("rm "+normfile);
	return lines;	
}

Array1dC<Point2dC> MBSplineC::LoadPointsFromFile(FilenameC const &f_name, UIntT const &size)
{
	IStreamC myfile(f_name);
	Array1dC<Point2dC> pts(size);
	UIntT index = 0;
	while((!myfile.eof())&&(index < size))
	{
		Point2dC pt(0.0,0.0);
		myfile>>pt.Row();
		myfile>>pt.Col();
		pts[index++] = pt.Copy();
	}
	myfile.Close();
	return pts;
}

bool MBSplineC::WritePointsToFile(Array1dC<Point2dC> const &pts, FilenameC const &f_name)
{
	ofstream myfile;
	myfile.open(f_name);
	if(myfile.is_open())
	{
		for(Array1dIterC<Point2dC> it(pts); it; it++)
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

void MBSplineC::CallSystem(StringC const &cmd)
{
	UIntT i = system(cmd);
	if(i != 0)
	{
		cout<<"Could not execute command - "<<cmd<<endl;
		exit(1);
	}
}
