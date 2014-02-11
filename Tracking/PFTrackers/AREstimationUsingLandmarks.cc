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
	//cout<<"Loaded Lip Data"<<endl;
	//MANUALLY LOAD ALL THE MODEL DATA
	//MEAN
	FilenameC mean_file = "/vol/vssp/lip-tracking/LipShapeModels/TrainingData/OmniMean.txt";
	IStreamC mean_str(mean_file);
	MatrixC mean(22,1);
	for(UIntT i = 0; i < 22; i++)
	{
		mean_str>>mean[i][0];
	}//cout<<"Mean"<<mean<<endl;
	//EVALS
	FilenameC evals_file = "/vol/vssp/lip-tracking/LipShapeModels/TrainingData/OmniEVals.txt";
	IStreamC evals_str(evals_file);
	MatrixC evals(22,1);
	for(UIntT i = 0; i < 22; i++)
	{
		evals_str>>evals[i][0];
	}//cout<<"Evals"<<evals<<endl;	
	//EVECT
	FilenameC evect_file = "/vol/vssp/lip-tracking/LipShapeModels/TrainingData/OmniEVect.txt";
	IStreamC evect_str(evect_file);
	MatrixC evect(22,22);
	for(UIntT i = 0; i < 22; i++)
	{
		for(UIntT j = 0; j < 22; j++)
		{
			evect_str>>evect[i][j];
		}
	}//cout<<"Evect"<<evect<<endl;	
	MatrixC evect_inv = evect.Inverse();	
	//cout<<"Completed PCA Models Loading"<<endl;
	//FOR AFFINE COMPUTATION
	RealT rows=142;	
  	RealT cols=120;
  	RealT rowFrac=0.35;
  	RealT colFrac=0.25;
  	Tuple4C<RealT,RealT,RealT,RealT> regParameters(rows,cols,rowFrac,colFrac); //For affine projection calculation
	DListC<StringC> qfile = qimg_dir.FiltList("*.ppm");
	//Perform eye_estimation for the first frame
	FilenameC firstimg = qimg_dir + qfile.First();
	Tuple2C<Point2dC,Point2dC> emergency_eyes = GetEyeLoc(firstimg);
	Array1dC<MatrixC> states(lp_coords.Size());
	//Array1dC<MatrixC> states(100);  
	UIntT ind = 0;
	for(DLIterC<StringC> it(qfile); it; it++)
	{
		if(ind < lp_coords.Size())
		//if(ind < 100)
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
			//Call Omni SW
			Tuple2C<Point2dC,Point2dC> eyes;
			try
			{
				eyes = GetEyeLoc(qimg);
			}
			catch(OmniN::ColossusExceptionNoFaceFoundC er)
			{
				eyes = emergency_eyes;
			}
			//Perform Affine Projection
			Tuple2C<ImageC<RealT>,Affine2dC> normface = convertImages(realimg, eyes.Data1(),eyes.Data2(), regParameters);
			ImageRectangleC outRect(Round(rows), Round(cols));
			WarpAffineC<RealT> warp(outRect,normface.Data2());
			ImageC<RealT> normalisedImage = warp.Apply(realimg);
			//cout<<"AFFINE WARP"<<endl;
			//B-Spline Stuff
			//Draw Red Crosses On Image
			//NOTE - ADDING AN EXTRA POINT - REPEATING FIRST POINT TO CREATE CLAMPED CURVE
			Array1dC<Point2dC> lp(lp_coords[ind].Size() + 1);
			for(UIntT i = 0; i < lp_coords[ind].Size(); i++)
			{
				lp[i].Row() = (lp_coords[ind])[i].Col();
				lp[i].Col() = (lp_coords[ind])[i].Row();
		   	//DrawCircle(currentimg,red,lp[i],5);    	
			}
			lp[lp_coords[ind].Size()] = lp[0];
			IntT order = 3;
			IntT ncp = 11;
			//Perform B-Spline Interpolation
			BSplineC bspl(order,ncp, BSplineC::UOPEN);
			Array1dC<Point2dC> cpts = bspl.CalculateControlPoints(lp, BSplineC::CHORDLENGTH, BSplineC::UOPEN,order, ncp);	
			Array1dC<Point2dC> ncpts = AffineControlPoints(cpts, normface.Data2().Inverse());
			//cout<<"NORMALISED INTERPOLATING CONTROL POINTS"<<endl;
			//Now we have our affine normalised control points, need to get eigenvalues
			VectorC lamda = GetEigenCoeff(ncpts,mean, evect_inv, numeig);
			//cout<<"EIGEN COEFFICIENT ESTIMATION"<<endl;
			VectorC aff_vect = GetAffAsVector(normface.Data2());//cout<<"AFFINE VECTOR"<<endl;
			//MatrixC state(VectorC(aff_vect.Append(lamda)));
			MatrixC state(lamda);
			states[ind] = state.Copy();
			ind++;
		}
	}
	//cout<<"States = "<<states<<endl;
	//NOW WE HAVE THE VALUES WE WANT, ESTIMATE THE AR COEFFICIENTS
	Array1dC<MatrixC> ar_coeff = Apply(states);
	cout<<"AR COEFFICIENTS \n"<<ar_coeff<<endl;
	return 0;
}

template <typename T>
DListC<Index2dC> ClusterIdent(DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > &finalpop, const ImageC<T> &img, const IntT &segtype)
{
	DListC<Index2dC> lip;
	//Method just chooses the type of function to use for cluster identification depending on segtype
	switch(segtype)
	{
		case 0:
		{
			lip = DensityCG(finalpop, img);
			break;
		}
		case 1:
		{
			lip = JCC(finalpop, img);
			break;
		}
	}
	return lip;
}


Tuple2C<ImageC<RealT>,Affine2dC> convertImages(const ImageC<RealT> & image, const Point2dC &nle,const Point2dC &nre, Tuple4C<RealT,RealT,RealT,RealT> regParameters)
{
    ImageC<RealT> normalisedImage;
/*
//!Gabor
    RealT rows = 160;
    RealT cols = 128;
//for leye=28c,55r
    RealT rowFrac = 0.34375;
    RealT colFrac = 0.21875;
*/
  RealT rows = regParameters.Data1();//142;
  RealT cols =  regParameters.Data2();//120;
  RealT rowFrac =  regParameters.Data3();//0.35;
  RealT colFrac =  regParameters.Data4();//0.25;

//for FRGC 3d Data
//        RealT rowFrac = 0.3;
 //   RealT colFrac = 0.1;
    UIntT gaussOrder = 7;
        
    //: Workout a simple geometric normalisation
    //==========================================
    Point2dC le(rows * rowFrac, (1.0-colFrac) * cols);
    Point2dC re(rows * rowFrac, cols * colFrac);
    ImageRectangleC outRect(Round(rows), Round(cols));
    Vector2dC d1 = le - re;
    Vector2dC d2 = nle - nre;
    RealT rot = d2.Angle() - d1.Angle();     
    RealT scale = d2.Modulus() / d1.Modulus();  
    Matrix2dC rotm = Matrix2dC(Cos(rot) * scale,-Sin(rot) * scale, Sin(rot) * scale,Cos(rot) * scale);      
    Point2dC cent= ((nle + nre)/2);
    Point2dC dcent = rotm * (le + re)/2;
    Vector2dC off = (cent - dcent);
    Affine2dC tr(rotm,off);
    
    //: Smooth and geometric normalise
    //=================================
    GaussConvolve2dC<RealT> smooth(gaussOrder);
    WarpAffineC<RealT> warp(outRect,tr);
    ImageRectangleC irec = warp.InputRectangle();
    irec = irec.Expand(gaussOrder/2 + 2);
    irec.ClipBy(image.Frame());
    normalisedImage = smooth.Apply(ImageC<RealT>(image,irec));

    normalisedImage = warp.Apply(normalisedImage); 
    //: Lets display the image
    //========================
	RavlN::Save("@X: normalised image", normalisedImage);
	Tuple2C<ImageC<RealT>,Affine2dC> out(normalisedImage, tr);
   return out;

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

Array1dC<Point2dC> AffineControlPoints(const Array1dC<Point2dC> &pts, const Affine2dC &aff)
{
	Array1dC<Point2dC> ncpts(pts.Size());
	//Apply affine transform to our points and draw onto normalised face image
	for(Array1dIter2C<Point2dC, Point2dC> it(pts,ncpts); it; it++)
	{
		Vector2dC vec(it.Data1().Row(), it.Data1().Col());
		vec = aff*vec;
		Point2dC ind(vec.Row(), vec.Col());
		it.Data2() = ind.Copy();  
	}
	return ncpts;
}

VectorC GetEigenCoeff(const Array1dC<Point2dC> &nc,const MatrixC &mu, const MatrixC &ev_inv, const IntT &num)
{
	//Convert the array into a vector
	VectorC vec(nc.Size()*2);
	vec.Fill(0);
	for(SizeT i = 0; i < nc.Size(); i++)
	{
		vec[i] = nc[i].Col();
		vec[i + nc.Size()] = nc[i].Row();
	}
	MatrixC cp_mat(vec);
	MatrixC evals = (cp_mat - mu)*ev_inv;
	//Convert matrix to vector
	VectorC res(num);
	for(IntT i = 0; i < num; i++)
	{
		res[i] = evals[i][0];
	}
	return res;
}

VectorC GetAffAsVector(const Affine2dC &af)
{
	VectorC res(6);
	res[0] = af.SRMatrix()[0][0];
	res[1] = af.SRMatrix()[0][1];
	res[2] = af.Translation()[0];
	res[3] = af.SRMatrix()[1][0];
	res[4] = af.SRMatrix()[1][1] ;
	res[5] = af.Translation()[1];
	return res;
}

MatrixC CalcRI(const IntT &ival, const Array1dC<MatrixC> &trainingdat)
{
	MatrixC out(trainingdat[0].Rows(),trainingdat[0].Cols());
	out.Fill(0.0);
	for(UIntT k = 2; k < trainingdat.Size(); k++)
	{
		out+= trainingdat[k-ival];
	}
	return out;
}

MatrixC CalcRIJ(const IntT &ival, const IntT &jval, const Array1dC<MatrixC> &trainingdat)
{
	MatrixC out(trainingdat[0].Rows(),trainingdat[0].Rows());
	out.Fill(0.0);
	for(UIntT k = 2; k < trainingdat.Size(); k++)
	{
		out+= (trainingdat[k-ival]).MulT(trainingdat[k-jval]);
	}
	return out;	
}

MatrixC CalcRCIJ(const IntT &ival, const IntT &jval, const Array1dC<MatrixC> &trainingdat)
{
	MatrixC out(trainingdat[0].Rows(),trainingdat[0].Rows());
	out.Fill(0.0);
	out = CalcRIJ(ival,jval,trainingdat);
	out -= (((CalcRI(ival,trainingdat).MulT(CalcRI(jval,trainingdat))))*(1.0/((RealT)trainingdat.Size() - 2.0)));
	return out;
}
 
Array1dC<MatrixC> Apply(const Array1dC<MatrixC> &tr)
{
	Array2dC<MatrixC> autocorr(3,3);
	for(IntT i = 0; i < 3; i++)
	{
		for(IntT j = 0; j < 3; j++)
		{
			autocorr[i][j] = CalcRCIJ(i,j,tr).Copy();
		}
	}
	///////////////////////
	MatrixC a2 = autocorr[0][2];
	MatrixC dummy = autocorr[1][1].Inverse();
	a2 -= ((autocorr[0][1]*dummy)*autocorr[1][2]);
	MatrixC a2sec = autocorr[2][2] - ((autocorr[2][1]*(MatrixC(autocorr[1][1].Inverse())))*autocorr[1][2]);
	a2 = a2*(a2sec.Inverse());
	MatrixC a1 = (((autocorr[0][1])-(a2)*(autocorr[2][1]))*((autocorr[1][1]).Inverse())).Copy();
	Array1dC<MatrixC> out(3);
	MatrixC d = (CalcRI(0,tr)-((a2)*(CalcRI(2,tr))) - ((a1)*(CalcRI(1,tr))));
	d = d*(1.0/((RealT)tr.Size() - 2.0));
	out[0] = a1.Copy();
	out[1] = a2.Copy();
	out[2] = d.Copy();
	return out;
}

