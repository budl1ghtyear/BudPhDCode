//ClosedBSplineFitting Program
#include 	"Ravl/SArray1d.hh"
#include 	"Ravl/SArray1dIter.hh"
#include 	"Ravl/SArray1dIter2.hh"
#include 	"Ravl/SArray1dIter3.hh"
#include 	"Ravl/Circle2d.hh" 
#include 	"Ravl/DArray1d.hh"
#include 	"Ravl/DArray1dIter.hh"
#include 	"Ravl/Image/DrawCross.hh"
#include 	"Ravl/Image/DrawLine.hh"
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
#include "BasisSpline.hh"
using namespace std;
using namespace RavlN;
using namespace RavlImageN;
using namespace RavlConstN;
//Simple Functions to Perform Some Data Generation
SArray1dC<Point2dC> GenerateCirclePoints(const UIntT &radius, const UIntT &npts); //Remember that the centre of the circle is at 50,50
void DisplayImage(const SArray1dC<Point2dC> &pts, const UIntT &size=2);
RealRGBValueC GenerateRandomColour(void);

RandomGaussC myrandgen;
int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	UIntT rad = opt.Int("r",100,"Circle Radius");
	UIntT npt = opt.Int("m",15,"Number of data Points to Generate");
	UIntT sze = opt.Int("s",5,"Size of display cross");
	//UIntT deg =  opt.Int("d",3,"Degree of the B-Spline");
	//UIntT npone = opt.Int("n",11,"Number of B-Spline control points (n+1)"); 
	opt.Check();
	SArray1dC<Point2dC> circ_pts = GenerateCirclePoints(rad,npt);
	cout<<"Generated Points are - "<<circ_pts<<endl;
	DisplayImage(circ_pts,sze);
	BasisSplineC spl(11,2,30);
	spl.FitCurve(circ_pts);
	cout<<"Fitted Control Points are - "<<spl.GetControlPoints()<<endl;
	cout<<"Fitted Curve Points are - "<<spl.GetCurvePoints()<<endl;
	cout<<"Fitted Curve Normals are - "<<spl.GetCurveNormals()<<endl;
	//sleep(5);
	spl.GenerateCurve(spl.GetControlPoints());
	cout<<"Generated Curve Points are - "<<spl.GetCurvePoints()<<endl;
	cout<<"Generated Curve Normals are - "<<spl.GetCurveNormals()<<endl;
	//Draw the resulting rendering:
	ImageC<RealRGBValueC> normalimg(IndexRange2dC(0,700,0,700));
	normalimg.Fill(RealRGBValueC(0.0,0.0,0.0));RealRGBValueC rgb = GenerateRandomColour();
	for(SArray1dIter2C<Point2dC, LinePP2dC> it(spl.GetCurvePoints(),spl.GetCurveNormals()); it; it++)
	{
		//DrawCross(normalimg,rgb,Index2dC(it.Data1().Row(),it.Data1().Col()),sze);
		Point2dC pt1 = it.Data2().Point(0.0); Index2dC from(pt1.Row(),pt1.Col());
		Point2dC pt2 = it.Data2().Point(10.0); Index2dC to(pt2.Row(),pt2.Col());
		cout<<"From - "<<from<<"\t To - "<<to<<"\t With vector = "<<it.Data2().Vector()<<endl;
		DrawLine(normalimg,rgb,from, to);     
	}
	if(!Save("@XA:Normal Image",normalimg)) cerr<<"Could not save the resulting output image"<<endl;
	return 0;
}

SArray1dC<Point2dC> GenerateCirclePoints(const UIntT &radius, const UIntT &npts)
{
	SArray1dC<Point2dC> pts(npts); 
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

void DisplayImage(const SArray1dC<Point2dC> &pts, const UIntT &size)
{
	ImageC<RealRGBValueC> img(IndexRange2dC(0,700,0,700));
	img.Fill(RealRGBValueC(0.0,0.0,0.0));
	//Generate Random Colour
	RealRGBValueC rgb = GenerateRandomColour();
	for(SArray1dIterC<Point2dC> it(pts); it; it++)
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
