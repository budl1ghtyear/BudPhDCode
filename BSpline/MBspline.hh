#ifndef MBSpline_HH
#define MBSpline_HH

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

/*
//HOME PATHS
StringC Mathematica = "/usr/local/bin/math";
StringC ScriptPath = "/home/bud/Desktop/MathematicaNB/";
StringC ResPath = "/home/bud/Desktop/SplineExperiments/";
*/

class MBSplineC
{
	public:
	//Default constructors
	MBSplineC(){}
	~MBSplineC(){}
	//Constructor
	MBSplineC(Array1dC<Point2dC> const &ctrlpts, UIntT const &deg=2 )
	{
		control_points = ctrlpts.Copy();
		degree = deg;
	}
	//Visible Methods
	Array1dC<Point2dC> FitQuadraticSpline(Array1dC<Point2dC> const &pts);
	Array1dC<Point2dC> RenderClosedQuadraticSpline(UIntT const &npts){return this->RenderClosedQuadraticSpline(this->control_points,npts);}
	Array1dC<Point2dC> RenderClosedQuadraticSpline(Array1dC<Point2dC> const &pts, UIntT const &npts);
	Array1dC<LinePP2dC> GetNormalLines(UIntT const &npts){return this->GetNormalLines(this->control_points,npts);}
	Array1dC<LinePP2dC> GetNormalLines(Array1dC<Point2dC> const &pts,  UIntT const &npts);
	LinePP2dC GetNormalLine(RealT const &par);
	//Accessor Methods
	Array1dC<Point2dC> GetControlPoints(void) const {return control_points;}
	void SetControlPoints(Array1dC<Point2dC> const &pts) {control_points = pts.Copy();}
	UIntT GetDegree(void) const {return degree;}
	void SetDegree(UIntT const &d) {degree = d;}
	static const StringC Mathematica;
	static const StringC ScriptPath;
	//StringC ResPath = "/vol/vssp/lip-tracking/SplineExperiments/";
	static const StringC ResPath;	
	protected:
	Array1dC<Point2dC> control_points;
	UIntT degree;
	Array1dC<Point2dC> LoadPointsFromFile(FilenameC const &f_name, UIntT const &size);
	bool WritePointsToFile(Array1dC<Point2dC> const &pts, FilenameC const &f_name);	
	void CallSystem(StringC const &cmd);
};

#endif
