//ConstructShapeModel Code
// Author - Bud Goswami
// Date - 9/11/09
// Process:
//1) Load the TRES DATA
//2) Perform either translation, rotation or scale normalisation on the TRES data
//3) Perform parameterisation of this raw data
//4) Output the data 

//Required Libs
#include "Ravl/Affine2d.hh"
#include "Ravl/Array2d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/SArray2d.hh"
#include "Ravl/SArray1dIter.hh"
#include "Ravl/DArray1d.hh"
#include "Ravl/DP/SequenceIO.hh"
#include "Ravl/DP/FileFormatIO.hh"
#include "Ravl/IO.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/DrawCircle.hh"
#include "Ravl/Image/DrawLine.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/LinePP2d.hh"
#include "Ravl/Math.hh"
#include "Ravl/Matrix.hh"
#include "Ravl/Option.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/SArray1d.hh"
#include "Ravl/SArray2d.hh"
#include "Ravl/SArray1dIter.hh"
#include "Ravl/Stream.hh"
#include "Ravl/Tuple2.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Vector2d.hh"
#include "Ravl/Vector.hh"
#include "Ravl/Image/WarpAffine.hh"
#include <iostream>
#include <fstream>
//#include <unistd>
#include "BSplineC.hh"
#include "ModellingUtilityFunctions.hh"
#include <cstdio>
#include "BasisSpline.hh"
using namespace RavlN;
using namespace RavlImageN;

Array1dC<Point2dC> SplineInterpolation(Array1dC<Point2dC> &data, Array1dC<RealT> &knots);
MatrixC ArrayToMatrix(const Array1dC<Point2dC> &pts);
Array1dC<Point2dC> MatrixToArray(const MatrixC &mat);
void DrawArray(ImageC<RealRGBValueC> &img, const Array1dC<Point2dC> &pts, const RealRGBValueC &rgb, const UIntT &size);
int main(int nargs,char **argv) 
{
	//CREATE THE REQUIRED INPUT OPTIONS
	OptionC opt(nargs,argv);
	//IO Options
	FilenameC lip = opt.String("l","in.txt","Text File with Lip-Co-ordinates");
	DirectoryC img_dir = opt.String("i","in.avi","Input stream. ");//specify input file
	IntT order = opt.Int("d",3,"Spline Order");
	IntT ncp = opt.Int("n",11,"Number of Control Points Required");
	opt.Check();
	
	ModellingUtilityFunctions mufc;
	//(1)Load the data
	DListC<SArray1dC<Point2dC> > lp_pts = mufc.LoadTRESData(lip);
	//Load the images
	DListC<StringC> lp_imgs = img_dir.FiltList("*.ppm");
	//Now iterate through both the lists and draw the resulting lip points as well as the B-spline control points
	DLIterC<SArray1dC<Point2dC> > pts_it(lp_pts);
	//Marker Colours:
	RealRGBValueC red(255,10,10),green(10,255,10),blue(10,10,255);
	for(DLIterC<StringC> img_it(lp_imgs); pts_it;img_it++)
	{
		//(2)Given some data, perform normalisation throughout
		Array1dC<Point2dC> unnormalised_raw_data = mufc.SArrayToArray(*pts_it);
		//DEBUG STEP: Draw the points
		FilenameC im_name = img_dir + (*img_it);
		ImageC<RealRGBValueC> img_draw((mufc.LoadImage(im_name)).Copy());
		DrawArray(img_draw,unnormalised_raw_data,red,2);
		
		//(3)Perform parameterisation of this raw data using B-Splines
		BSplineC bspl(order,ncp, BSplineC::UPERIODIC);
		Array1dC<Point2dC> lip_cp = bspl.CalculateControlPoints(unnormalised_raw_data,BSplineC::UOPEN,BSplineC::UNIFORM);
		cout<<"Control Points for BSPline = "<<lip_cp<<endl;
		DrawArray(img_draw,lip_cp,green,4);
		Array1dC<Point2dC> spl_pts = bspl.RenderCurve(300);
		DrawArray(img_draw,spl_pts,blue,1);
		BasisSplineC bas(order,ncp);
		Array1dC<RealT> openuniformknotvector = bas.OpenUniformKnotVector();
		//cout<<"Open Uniform Knot Vector "<<openuniformknotvector<<"\n Par_range = "<<bas.GetParameterRange().Data1()<<"\t"<<bas.GetParameterRange().Data2()<<endl;
		Array1dC<RealT> parameterisation = bas.Parameterisation(unnormalised_raw_data,BasisSplineC::UNIFORM);
		//cout<<"Parameters - "<<parameterisation;
		//~ Array1dC<RealT> periodicuniformknotvector = bas.PeriodicUniformKnotVector();
		//~ cout<<"Periodic Uniform Knot Vector "<<periodicuniformknotvector<<"\n Par_range = "<<bas.GetParameterRange()<<endl;
		Array1dC<RealT> parameterisation2 = bas.Parameterisation(unnormalised_raw_data,BasisSplineC::CHORDLENGTH);
		cout<<"Chord Length Parameters - "<<parameterisation2;
		cout<<"Order = "<<order<<"\t Number of Control Points = "<<ncp<<endl;
		//Array1dC<RealT> basis = bas.ComputeBasisFunction(order, parameterisation[1], ncp,openuniformknotvector);
		//cout<<"Basis = "<<basis<<endl;
		Array1dC<Point2dC> control_poly = bas.BSplineFit(unnormalised_raw_data,BasisSplineC::UNIFORM, BasisSplineC::OPENUNIFORM);
		cout<<"Control Polygon = "<<control_poly<<endl;
		cout<<"Basis matrix = "<<bas.ComputeNMatrix()<<endl;
		//Closed Averaging Example:
		Array1dC<RealT> knotavegclosed = bas.ClosedAveragedKnotVector(parameterisation2);
		cout<<"Closed Averaged Knot Vector = "<<knotavegclosed<<endl;
		Array1dC<IntT> knotspans(parameterisation2.Size()); knotspans.Fill(0);
		for(Array1dIter2C<RealT, IntT> it(parameterisation2,knotspans); it; it++)
		{
			it.Data2() = bas.FindSpan(parameterisation2,knotavegclosed,it.Data1());
		}
		cout<<"Knot spans - "<<knotspans<<endl;
		DrawArray(img_draw,control_poly,blue,5);
		sleep(100);
	}
	return 0;	
}

void DrawArray(ImageC<RealRGBValueC> &img, const Array1dC<Point2dC> &pts, const RealRGBValueC &rgb, const UIntT &size)
{
	for(Array1dIterC<Point2dC> it(pts); it; it++)
	{
		DrawCross(img,rgb,(*it),size);
	}
	if(!Save("@XA:Drawn Image",img)) cerr<<"Could not show the drawing"<<endl;
}
//Perform B-Spline interpolation according to our book
Array1dC<Point2dC> SplineInterpolation(Array1dC<Point2dC> &data, Array1dC<RealT> &knots)
{
	MatrixC dat = ArrayToMatrix(data);
	RealT t1=0.0,t2=0.0,t3=0.0,t4=0.0,t5=0.0;
	UIntT nodata = data.Size(); cout<<"Num data = "<<nodata<<endl;
	Array1dC<RealT> alpha(nodata),beta(nodata), gamma(nodata);
	alpha[0] = 0; beta[0] = 1; gamma[0] = 0;
	alpha[nodata - 1] = 0; beta[nodata - 1] = 1; gamma[nodata - 1] = 0; 
	for(UIntT i = 0; i < nodata -1; i++)
	{
			t1 = knots[i+1];
			t2 = knots[i+2];
			t3 = knots[i+3];
			t4 = knots[i+4];
			t5 = knots[i+5];
			alpha[i] = (t4 - t3)*(t4 - t3)/(t4 - t1);
			beta[i] = (t3 - t1)*(t4 - t3)/(t4 - t1) + (t5 - t3)*(t3 - t2)/(t5 - t2);
			gamma[i] = (t3 - t2)*(t3 - t2)/(t5 - t2);
			if((t4 - t2)== 0.0) //To prevent division by zero
			{
				alpha[i] = 0.0;
				beta[i] = 0.0;
				gamma[i]= 0.0 ;
			}
			else
			{
				alpha[i] /= (t4 - t2);
				beta[i] /= (t4 - t2);
				gamma[i] /= (t4 - t2);
			}
	}
	cout<<"Alpha = "<<alpha;cout<<"Beta = "<<alpha;cout<<"Gamma = "<<alpha;
	return data;
}

MatrixC ArrayToMatrix(const Array1dC<Point2dC> &pts)
{
	MatrixC mat(pts.Size(),2);
	for(UIntT i = 0; i < pts.Size(); i++)
	{
		mat[i][0] = pts[i].Row();
		mat[i][1] = pts[i].Col();
	}
	return mat;
}

Array1dC<Point2dC> MatrixToArray(const MatrixC &mat)
{
	Array1dC<Point2dC> arr(mat.Rows());
	for(UIntT i = 0; i < mat.Rows(); i++)
	{
		Point2dC pt(mat[i][0],mat[i][1]);
		arr[i] = pt.Copy();
	}
	return arr;
}

