//ClosedBSplineFitting Program
#include 	"Ravl/Array1d.hh"
#include 	"Ravl/Array1dIter.hh"
#include 	"Ravl/Array1dIter3.hh"
#include 	"Ravl/Circle2d.hh" 
#include 	"Ravl/DArray1d.hh"
#include 	"Ravl/DArray1dIter.hh"
#include 	"Ravl/Image/DrawCross.hh"
#include 	"Ravl/Image/Image.hh"
#include 	"Ravl/Image/RealRGBValue.hh"
#include 	"Ravl/IndexRange2d.hh"
#include 	"Ravl/IO.hh"
#include 	"Ravl/LinePP2d.hh"
#include 	"Ravl/Option.hh"
#include 	"Ravl/Point2d.hh"
#include 	"Ravl/RandomGauss.hh"
#include 	"Ravl/RandomMersenneTwister.hh"
#include 	"Ravl/Vector2d.hh"
#include 	<iostream>
#include 	<fstream>
#include 	"Ravl/OS/Filename.hh"
#include 	<stdio.h>
#include 	<stdlib.h>
#include 	"Ravl/String.hh"

using namespace std;
using namespace RavlN;
using namespace RavlImageN;
using namespace RavlConstN;
//Simple Functions to Perform Some Data Generation
Array1dC<Point2dC> GenerateCirclePoints(const UIntT &radius, const UIntT &npts); //Remember that the centre of the circle is at 50,50
void DisplayImage(const Array1dC<Point2dC> &pts, const UIntT &size=2);
RealRGBValueC GenerateRandomColour(void);

Array1dC<Point2dC> FitQuadraticSpline(Array1dC<Point2dC> const &pts);
Array1dC<Point2dC> RenderClosedQuadraticSpline(Array1dC<Point2dC> const &pts, UIntT const &npts);
Array1dC<LinePP2dC> GetNormalLines(Array1dC<Point2dC> const &pts, UIntT const &npts);
Array1dC<Point2dC> LoadPointsFromFile(FilenameC const &f_name, UIntT const &size);
bool WritePointsToFile(Array1dC<Point2dC> const &pts, FilenameC const &f_name);

RandomGaussC myrandgen;
//~ StringC Mathematica = "/usr/local/bin/math";
//~ StringC ScriptPath = "/home/bud/Desktop/MathematicaNB/";
//~ StringC ResPath = "/home/bud/Desktop/SplineExperiments/";
//WORK PATHS
StringC Mathematica = "/vol/vssp/lip-tracking/ExtSW/Mathematica/Scripts/math";
StringC ScriptPath = "/vol/vssp/lip-tracking/MathematicaNB/";
StringC ResPath = "/vol/vssp/lip-tracking/SplineExperiments/";

int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	UIntT rad = opt.Int("r",100,"Circle Radius");
	UIntT npt = opt.Int("m",15,"Number of data Points to Generate");
	UIntT sze = opt.Int("s",5,"Size of display cross");
	//UIntT deg =  opt.Int("d",3,"Degree of the B-Spline");
	//UIntT npone = opt.Int("n",11,"Number of B-Spline control points (n+1)"); 
	opt.Check();
	Array1dC<Point2dC> circ_pts = GenerateCirclePoints(rad,npt);
	cout<<"Generated Points are - "<<circ_pts<<endl;
	DisplayImage(circ_pts,sze);
	Array1dC<Point2dC> c_pts = FitQuadraticSpline(circ_pts);
	DisplayImage(c_pts,3);
	Array1dC<Point2dC> curve_pts = RenderClosedQuadraticSpline(c_pts, 30); cout<<"curve_pts "<<curve_pts.Size()<<endl;
	DisplayImage(curve_pts,10);
	Array1dC<LinePP2dC> normlines = GetNormalLines(c_pts,30);
	
	return 0;
}

Array1dC<Point2dC> GenerateCirclePoints(const UIntT &radius, const UIntT &npts)
{
	Array1dC<Point2dC> pts(npts); 
	Point2dC center(350,350);
	Circle2dC circle(center,radius);
	RealT incr = (2.0*pi)/(RealT)npts;
	for(UIntT i = 0; i < npts; i++)
	{
		RealT angle = (RealT)i * incr;
		Point2dC pt = circle.Value(angle);
		pts[i] = pt.Copy();
	}
	return pts;
}

void DisplayImage(const Array1dC<Point2dC> &pts, const UIntT &size)
{
	ImageC<RealRGBValueC> img(IndexRange2dC(0,700,0,700));
	img.Fill(RealRGBValueC(0.0,0.0,0.0));
	//Generate Random Colour
	RealRGBValueC rgb = GenerateRandomColour();
	for(Array1dIterC<Point2dC> it(pts); it; it++)
	{
		DrawCross(img,rgb,Index2dC((*it).Row(),(*it).Col()),size);
	}
	if(!Save("@XA:Output Image",img)) cerr<<"Could not save the resulting output image"<<endl;
}

RealRGBValueC GenerateRandomColour(void)
{
	//Generate 3 Random Numbers Between 1 and 255
	RealRGBValueC col((RealT)(myrandgen.Generate()*255),(RealT)(myrandgen.Generate()*255),(RealT)(myrandgen.Generate()*255));
	return col;
}

Array1dC<Point2dC> FitQuadraticSpline(Array1dC<Point2dC> const &pts)
{
	//Write the points to an output file just as space separated points x1 y1 \n x2 y2
	FilenameC ptfile = ResPath + "000LipPoints.txt";
	WritePointsToFile(pts,ptfile);
	//Create a command file
	StringC cmd = Mathematica + " -run \"<< " + ScriptPath + "make_control_points.m\"";
	system(cmd);
	//Now load in the input file
	FilenameC resfile = ResPath + "000LipPoints.txt.csv";
	Array1dC<Point2dC> ctrlPts = LoadPointsFromFile(resfile,11);
	//Clean up residuals
	StringC rmcmd1 = "rm "+ptfile; system(rmcmd1);
	StringC rmcmd2 = "rm "+resfile; system(rmcmd2);	
	return ctrlPts;
}

Array1dC<Point2dC> RenderClosedQuadraticSpline(Array1dC<Point2dC> const &pts, UIntT const &npts)
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
	system(cmd);
	//Now load in the resulting input file
	FilenameC resfile = ResPath + "000CurvePoints.txt";
	Array1dC<Point2dC> curvePts = LoadPointsFromFile(resfile,npts);
	//Clean up residuals
	//system("rm "+resfile);system("rm "+ptfile);system("rm " + numpt);
	return curvePts;
}

Array1dC<LinePP2dC> GetNormalLines(Array1dC<Point2dC> const &pts, UIntT const &npts)
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
	system(cmd);
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
	system("rm "+resfile);system("rm "+ptfile);system("rm " + numpt);
	return lines;	
}

Array1dC<Point2dC> LoadPointsFromFile(FilenameC const &f_name, UIntT const &size)
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

bool WritePointsToFile(Array1dC<Point2dC> const &pts, FilenameC const &f_name)
{
	ofstream myfile;
	myfile.open (f_name);
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

