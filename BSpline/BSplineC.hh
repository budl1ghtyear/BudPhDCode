//File - BSplineC.hh
// Author - Bud Goswami
// Date - 1/04/09
// Function - to define a B-Spline Curve class
// Please Note that the functionality in this class including variable naming conventions are adapted from:
// CH3 of the book - Introduction to NURBS: A Historical Perspective by Dave Rogers

#ifndef BSplineC_H
#define BSplineC_H
//////////////////////////////////////////////////////////////////////////////
//! file =      "Ravl/LipCode/BSpline/BSplineC.hh"   
//! author =    "Bud Goswami"
//! lib =       SegLib
//! date =      "01/04/09"
//! userlevel = Normal
//! docentry =  "Ravl.API.LipCode.BSpline"
//Include Files
#include 	"Ravl/Array1d.hh"
#include 	"Ravl/Array1dIter.hh"
#include 	"Ravl/Array1dIter2.hh"
#include 	"Ravl/Index2d.hh"
#include  	"Ravl/Point2d.hh"
#include 	"Ravl/Math.hh"
#include 	"Ravl/Matrix.hh"
#include 	"Ravl/Point2d.hh"
#include 	"Ravl/LinePP2d.hh"
#include 	"Ravl/LineABC2d.hh"
#include 	"Ravl/Vector2d.hh"
#include 	"Ravl/Math.hh"
//Required Libraries
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace RavlN {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------------------------------------------------------------------
// **********  BSplineC  **********************************************
// ---------------------------------------------------------------------------
// (previous lines will be ignored - no preceeding //: )

//: This is a BSpline Curve Class, written with the aim of performing curve fitting, curve rendering and normal drawing for a set of data co-ordinates
//
// Please Note that the functionality in this class including variable naming conventions are adapted from:<p>
// CH3 of the book - Introduction to NURBS: A Historical Perspective by Dave Rogers<p>

class BSplineC
{
public:
	//ENUM constant can be applied here
	// Required to enable the user to define parameter selection and knot vector generation
	enum KVG {UPERIODIC,UOPEN,NPERIODIC, NOPEN,UPCLOSED}; 
	//:Knot vector generation - Uniform Periodic, Uniform Open, Non-uniform Peridic and Non-uniform Open
	enum ParSel{UNIFORM,CHORDLENGTH,CENTRIPETAL}; 
	//:Parameter Selection - Uniformly space, chord-length, centripetal

	//Constructors:
	BSplineC(){}
	~BSplineC(){}
	BSplineC(const UIntT &ord, const Array1dC<RealT> &kv,Array1dC<Point2dC> cp);
	//: Constructor that takes in the order, knot vector and control point vector
	//BSplineC(const UIntT &ord, const UIntT &ncp, KVG ktype);
	BSplineC(const UIntT &ord, const UIntT &ncp, KVG ktype, bool status=true);
	//: Constructor that takes in the order, number of control points, Knot vector generation type
	//Obtain Parameter Vectors
	Array1dC<RealT> GetParameters(const Array1dC<Point2dC> &points, ParSel ptype);
	//: Perform parameter generation using some data points and a specified method of parameter selection
	//Obtain Knot Vector
	Array1dC<RealT> CalculateKnotVector(const UIntT &ord, const UIntT &nc, KVG ktype);
	//: Calculate the knot vector using the order, number of control points and knot vector generation type
	Array1dC<RealT> CalculateKnotVector(const UIntT &ord, const UIntT &nc, const Array1dC<Point2dC> &cp, KVG ktype);
	//: Calculate the knot vector using the order, number of control points, control point vector and knot vector generation type
	Array1dC<Point2dC> CalculateControlPoints(const Array1dC<Point2dC> &points, KVG ktype, ParSel ptype);
	//: Calculate the control points using some data points, a knot vector generation type and a parameter selection type
	Array1dC<Point2dC> CalculateControlPoints(const Array1dC<Point2dC> &points, KVG ktype, ParSel ptype, const UIntT &ord, const UIntT &ncp);
	//: Calculate the control points using some data points, a knot vector generation type, a parameter selection type, B-Spline order and number of control points
	Array1dC<Point2dC> CalculateControlPoints(const Array1dC<Point2dC> &points, const Array1dC<RealT> &kv, ParSel ptype);
	//: Calculate the control points using some data points, a knot vector, a parameter selection type
	
	//Accessor Method for KV
	Array1dC<RealT> GetKnotVector(void) const
	//: Accessor method that returns the computed knot vector
	{	
		//if(this->knot_vec.Size() == 0)
		//{
		//	CalculateKnotVector(this->order, this->numcp, UOPEN);
		//	cout<<"knot_vec has as yet not been computed. Computing UNIFORM OPEN knot_vector"<<endl;
		//}
		return knot_vec;
	}
	void SetKnotVector(const Array1dC<RealT> &kv) {knot_vec = kv.Copy();}
	//: Set knot vector with generated values
	Array1dC<Point2dC> CalculateControlPoints(const Array1dC<Point2dC> &points, ParSel ptype, KVG ktype,const UIntT &ord, const UIntT &ncp);
	
	//Accessor Method for Cp
	Array1dC<Point2dC> GetControlPoints(void) const {return ctrl_vec;}
	//: Accessor method for control point vector
	void SetControlPoints(const Array1dC<Point2dC> &pt) 
	//: Accessor method to set the control point vector
	{
		ctrl_vec = pt.Copy();
	}
	//Other Accessor Methods:
	UIntT GetOrder(void) const {return order;}
	//: Get the order of the B-Spline
	void SetOrder(const UIntT &ord)
	//: Set the order of the B-Spline
	{
		if((ord >= 2) && (numcp == 0)) //check if number of control points has been set
		{
			order = ord;
		}
		else if((ord >= 2) && (numcp != 0)) //if number of control points has been specified
		{
			if((ord >= 2) && (ord <= numcp))
				order = ord;
			else
				cerr<<"ERROR - ORDER, k must satisfy condition: 2 <= k <= (n+1), where (n+1) = number of Control Points Desired"<<endl;
		}	
		else
			cerr<<"ERROR - ORDER, k must satisfy condition: 2 <= k <= (n+1), where (n+1) = number of Control Points Desired"<<endl;
	} 
	UIntT GetNumControlPoints(void) const {return numcp;}
	//: Get number of control points
	void SetNumControlPoints(const UIntT &nc)
	//: Set number of control points
	{
		if(nc > 2)
			numcp = nc;
		else
			cerr<<"ERROR - Number of control poinst MUST be greater than 2 for a B-Spline";
	}	
	UIntT GetDegree(void) const{ return degree;} 
	//: Get the degree of the B-Spline. No SetDegree Function (degree = order - 1)
	Array1dC<Point2dC> RenderCurve(const UIntT &res);
	//: Returns the points that are generated using a curve rendering algorithm using a resolution of "res". This means that the parameter value will range from 0 to 1 with increments of 1/res
	Array1dC<Point2dC> RenderClosedCurve(const UIntT &res);
	//: Render closed curve. This method may not be functioning as expected
	LinePP2dC GetNormalLine(const RealT &par) const;
	//: Get the line equation specified by the normal to the curve at some parameter value (this val ranges from 0 to 1)
	Point2dC GetPointOnCurve(const RealT &par) const;
	//: Get a specific point on the actual B-Spline curve specified by some parameter value
protected:
	//Data Variables
	bool open;
	UIntT order; 
	//:K
	UIntT degree; 
	//:K-1
	UIntT numcp; 
	//:N+1
	Array1dC<RealT> knot_vec; 
	//:Knot Vector
	Array1dC<Point2dC> ctrl_vec; 
	//:Control Points
	//Parameter Selection Methods:
	Array1dC<RealT> ChordLengthParSel(const Array1dC<Point2dC> &points);
	//: Perform chorg length parameter selection
	Array1dC<RealT> CentripetalParSel(const Array1dC<Point2dC> &points);
	//: Perform centripetal method of parameter selection
	Array1dC<RealT> UniformParSel(const Array1dC<Point2dC> &points);
	//Perform uniform method of parameter selection
	//Knot Vector Generation Methods
	Array1dC<RealT> UniformPeriodicKVG(const UIntT &ord, const UIntT &nc);
	//: Uniform periodic knot vector
	Array1dC<RealT> UniformOpenKVG(const UIntT &ord, const UIntT &nc);
	//: Uniform open knot vector
	Array1dC<RealT> NonUniformPeriodicKVG(const UIntT &ord, const UIntT &nc, const Array1dC<Point2dC> &cp);
	//: Non-uniform periodic knot vector
	Array1dC<RealT> NonUniformOpenKVG(const UIntT &ord, const UIntT &nc, const Array1dC<Point2dC> &cp);	
	//: Non-uniform open knot vector
	Array1dC<RealT> ClosedUniformPeriodicKVG(const UIntT &ord, const UIntT &nc);	
	//: Closed uniform periodic knot vector. Closed curve gen at the moment does not function as expected. 
	Array1dC<RealT> NormaliseKnotVector(const Array1dC<RealT> &kv);
	//: Normalise the knot vector so it ranges from [0 to 1]
	//Basis Function generator
	Array1dC<RealT> ComputeBasis(const UIntT &ord, const UIntT &nc, const RealT &par, const Array1dC<RealT> &kv) const;
	//: Compute the deBoor basis using the deBoor recursion
	MatrixC GetBasisMatrix(const UIntT &ord, const UIntT &nc, const Array1dC<RealT> &pars, const Array1dC<RealT> &kv);
	//: For some parameter values, compute a basis matrix
	Array1dC<RealT> ComputeFirstDer(const UIntT &ord, const UIntT &nc, const RealT &par, const Array1dC<RealT> &kv) const;	
	//:Compute First Derivative
	MatrixC GetMatrixFromArray(const Array1dC<Point2dC> &points);
	//:Simple methods for data structure conversion
};
}
#endif
