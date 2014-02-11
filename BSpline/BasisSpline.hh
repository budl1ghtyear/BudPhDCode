#ifndef BasisSpline_HH
#define BasisSpline_HH


//This class is a wrapper class to a mathematica interface that performs closed B-Spline interpolation and related curve drawing operations
#include 	"Ravl/SArray1d.hh"
#include 	"Ravl/SArray1dIter.hh"
#include 	"Ravl/SArray1dIter3.hh"
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
//WORK PATHS

class BasisSplineC
{
	public:
	//Default constructors
	
	BasisSplineC(UIntT const &ncp=11, UIntT const &deg = 2, UIntT const &num = 30);
	//Constructor to use when perform curve generation
	BasisSplineC(SArray1dC<Point2dC> const &ctrlpts, UIntT const &deg = 2, UIntT const &num = 30);
	BasisSplineC(SArray1dC<Point2dC> const &cpts,SArray1dC<Point2dC> const &curvpts,SArray1dC<LinePP2dC> const &cnrm, UIntT const &ncp = 11, UIntT const &deg = 2, UIntT const &npts = 30 ):controlpoints(cpts),curvepoints(curvpts),curvenormals(cnrm),numcp(ncp),degree(deg),numpts(npts) {;}
	BasisSplineC Copy() const
	{
		return BasisSplineC(controlpoints,curvepoints,curvenormals,numcp, degree, numpts);
	}
	//Accessor Functions: SHOULD ONLY BE CALLED AFTER AN INSTANCE OF APPLY HAS BEEN CALLED
	SArray1dC<Point2dC> GetControlPoints(void) const {return controlpoints;}
	SArray1dC<Point2dC> GetCurvePoints(void) const {return curvepoints;}
	SArray1dC<LinePP2dC> GetCurveNormals(void) const {return curvenormals;}
	
	void FitCurve(SArray1dC<Point2dC> const &pts);
	void GenerateCurve(SArray1dC<Point2dC> const &ctrlpts);	
	static const StringC Mathematica;
	static const StringC ScriptPath;
	//StringC ResPath = "/vol/vssp/lip-tracking/SplineExperiments/";
	static const StringC ResPath;
		
	protected:
	BasisSplineC(){}
	SArray1dC<Point2dC> controlpoints; //control points
	SArray1dC<Point2dC> curvepoints; //rendered curve points
	SArray1dC<LinePP2dC> curvenormals; //the collection of generated curve normals
	UIntT numcp; //number of control points
	UIntT degree; //degree of the interpolation B-Spline
	UIntT numpts; //number of points for curve rendering and the curve values at those points
	
	bool LoadPointsFromFile(FilenameC const &f_name, UIntT const &size, SArray1dC<Point2dC> &res);
	bool WritePointsToFile(SArray1dC<Point2dC> const &pts, FilenameC const &f_name);	
	bool WriteIntToFile(UIntT const &n, FilenameC const &f_name);
	void CallSystem(StringC const &cmd);
	void RemoveFiles(FilenameC const &file);
};

#endif
