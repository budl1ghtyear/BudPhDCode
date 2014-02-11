//OmniSOAREstimate
//Author - Bud Goswami
//Date - 20 April 09
//Purpose: to perform 2nd Order AR estimation using a training video
//Process:
//Inputs - training video, mean, eigenvalues, eigenvectors, landmarked lip data
//Perform Lip Segmentation Using CMCD, NormRG
//Extract the B-Spline Lip Contour and save into a matrix
//For each step, convert this B-Spline matrix with eye locations into a statevector that is saved
//Perform AR estimation for that subject using these estimates and return onto screen

//INCLUDE FILES
#include "Ravl/Array1d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/Array2dIter2.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/DP/SequenceIO.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/DrawLine.hh"
#include "Ravl/Ellipse2d.hh"
#include "Ravl/Eigen.hh"
#include "Ravl/Image/DeinterlaceStream.hh"
#include "Ravl/Image/DrawCircle.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/DrawLine.hh"
#include "Ravl/Image/EdgeSobel.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Index2d.hh"
#include "Ravl/IndexRange2d.hh"
#include "Ravl/IO.hh"
#include "Ravl/LinePP2d.hh"
#include "Ravl/MeanCovariance.hh"
#include "Ravl/Option.hh"
#include "Ravl/PatternRec/DesignFuncPCA.hh"
#include "Ravl/PatternRec/Function.hh"
#include "Ravl/PatternRec/FuncMeanProjection.hh"
#include "Ravl/PatternRec/Sample.hh"
#include "Ravl/Random.hh"
#include "Ravl/StateVector.hh"
#include "Ravl/Tuple2.hh" 
#include "Ravl/VectorMatrix.hh"
//USER DEFINED CLASSES
#include "BSplineC.hh"
#include "ColourConvert.hh"//Contains colour conversion functions
#include "MouthRegion.hh"	//Contains mouth-region identification Functions 
#include "CostFunctions.hh"//Contains functions related to cluster grouping
#include "CascadedMCD.hh"
#include "BasePropagationModel.hh"
RandomGaussC BasePropagationModelC::rnd;
Tuple2C<ImageC<RealT>,Affine2dC> convertImages(const ImageC<RealT> & image, const Point2dC &nle,const Point2dC &nre, Tuple4C<RealT,RealT,RealT,RealT> regParameters);
ImageC<RealT> ConvertRealRGBToReal(const ImageC<RealRGBValueC> &rgb);
//Compute Affine Transform Projection on Control Points
Array1dC<Point2dC> AffineControlPoints(const Array1dC<Point2dC> &pts, const Affine2dC &aff);
VectorC GetAffAsVector(const Affine2dC &af);
//Compute EigenCoefficients of the data
VectorC GetEigenCoeff(const Array1dC<Point2dC> &nc,const MatrixC &mu, const MatrixC &ev_inv, const IntT &num);
MatrixC CalcRI(const IntT &ival, const Array1dC<MatrixC> &trainingdat);
MatrixC CalcRIJ(const IntT &ival, const IntT &jval, const Array1dC<MatrixC> &trainingdat);
MatrixC CalcRCIJ(const IntT &ival, const IntT &jval, const Array1dC<MatrixC> &trainingdat);
Array1dC<MatrixC> Apply(const Array1dC<MatrixC> &tr);
int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	StringC lip = opt.String("l","in.txt","Text File with Lip-Co-ordinates");
	DirectoryC qimg_dir = opt.String("i","TestSet/","Input Query Image Directory"); //For visualisation process
	FilenameC cascade_name = opt.String("c","/usr/share/opencv/haarcascades/haarcascade_frontalface_alt.xml","Path for the Haar Classifier");
	IntT numeig = opt.Int("g",5,"Number of EigenCoefficients To Use");
	opt.Check();
	
	//LOAD LIP DATA
	DArray1dC<SArray1dC<Point2dC> > lp_coords;
	IStreamC is(lip);
	while(!is.IsEndOfStream())
	{
		SArray1dC<Point2dC> pts;
		is >> pts;
		if(pts.Size() > 0)
		{
			lp_coords.Append(pts);			
		}		
	}
	DListC<StringC> qfile = qimg_dir.FiltList("*.ppm");
	//Perform eye_estimation for the first frame
	FilenameC firstimg = qimg_dir + qfile.First();
	Tuple2C<Point2dC,Point2dC> emergency_eyes = GetEyeLoc(firstimg);
	UIntT ind = 0;
	for(DLIterC<StringC> it(qfile); it; it++)
	{
		if(ind < lp_coords.Size())
		{
			FilenameC qimg = qimg_dir+(*it);
			ImageC<RealRGBValueC> src;   
			if(!Load(qimg, src)) cerr<<"Loading RAVL Image Failed"<<endl;
			FilenameC inp_file = qimg_dir + (*it);
			//Get Current Image
			ImageC<RealRGBValueC> currentimg;
			if(!Load(inp_file, currentimg)) cerr<<"Input file could not be loaded"<<endl;
			//Convert IMAGE to Real
			ImageC<RealT> realimg = ConvertRealRGBToReal(currentimg);
			Array1dC<Point2dC> lps = lp_coords[ind++];
			for(UIntT i = 0; i < lps.Size()-1; i++)
			{
				DrawLine(realimg,0.0,lps[i],lps[i+1]);   
			}
			if(!Save("@X:Lip Region",realimg)) cerr<<"Could not output the lip image"<<endl;
		}
	}
	return 0;
}

ImageC<RealT> ConvertRealRGBToReal(const ImageC<RealRGBValueC> &rgb)
{
	ImageC<RealT> img(rgb.Frame(),0.0);
	for(Array2dIter2C<RealRGBValueC, RealT> it(rgb,img); it; it++)
	{
		it.Data2() = it.Data1().Y();
	}	
	return img;
}


