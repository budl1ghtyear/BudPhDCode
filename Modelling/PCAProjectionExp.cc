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
#include "Ravl/PatternRec/FuncMeanProjection.hh"
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
#include "MBspline.hh"
#include "ModellingUtilityFunctions.hh"
#include "Align.hh"
#include <cstdio>
#include "ASMRotationNormalisation.hh"
#include "ASMAffineNormalisation.hh"
#include "Ravl/PatternRec/SampleIter.hh"
#include "Ravl/Image/AAMAppearance.hh"
#include "Ravl/RandomGauss.hh"
#include "Ravl/RandomMersenneTwister.hh"
#include "PCAProjection.hh"
#include "SimilarityProjection.hh"
#include "ICAProjection.hh"
using namespace RavlN;
using namespace RavlImageN;

RandomGaussC myrandgen;
VectorC PerformPCAProjection(VectorC const &mu, VectorC const &lamda, MatrixC const &e);
SampleC<SArray1dC<Point2dC> > PerformComponentProjection(VectorC const &mu, VectorC const &vals, MatrixC const &e, UIntT const &comp=0, IntT const &dev=0, UIntT const &ncp = 22 );
SArray1dC<Point2dC> VectorToSArray(VectorC const &vec);
void DisplayImage(const SArray1dC<Point2dC> &pts, ImageC<RealRGBValueC> &img, const UIntT &size=2);
RealRGBValueC GenerateRandomColour(void);
int main(int nargs,char **argv) 
{
	//CREATE THE REQUIRED INPUT OPTIONS
	OptionC opt(nargs,argv);
	//IO Options
	FilenameC lip = opt.String("d","in.txt","Text File with Lip-Co-ordinates in SampleC<SArray1dC<Point2dC> > stream format");
	DirectoryC pca = opt.String("p","./PCA/","Directory containing the PCA data");
	DirectoryC icadir = opt.String("i","./JADE/","Directory containing the ICA data");
	StringC ext = opt.String("e","spl","Extension of the file containing the lip-coordinates in SampleC<SArray1dC<Point2dC> > format");
	UIntT eigencomponent = opt.Int("c",1,"EigenComponent Along Which to Perform Projection");
	IntT devs = opt.Int("s",1,"Number of standard deviations to project");
	UIntT numcomponents = opt.Int("n",10,"Number of principal components to keep");
	opt.Check();
	/*
	//Read in the necessary PCA functions
	FilenameC pmean = pca+"Mean.tsv";
	VectorC pca_mean;
	if(LoadFileFromStream(pmean,pca_mean))
		cout<<"Loaded the pca mean - "<<pca_mean<<endl;
	
	FilenameC psdev = pca+"Eigenvalues.tsv";
	VectorC pca_sdev;
	if(LoadFileFromStream(psdev,pca_sdev))
		cout<<"Loaded the pca eigenvalues - "<<pca_sdev<<endl;	
	
	FilenameC ploadings = pca+"Eigenvectors.tsv";
	MatrixC pca_loadings;
	if(LoadFileFromStream(ploadings,pca_loadings))
		cout<<"Loaded the pca eigenvectors - "<<pca_loadings<<endl;	
	
	FilenameC affinemean = pca+"MeanPoints.txt";
	SArray1dC<Point2dC> aff_mean;
	if(LoadFileFromStream(affinemean, aff_mean))
		cerr<<"Affine Mean - "<<aff_mean<<endl;	
	
	//Now simply perform sdev projection using two parameters - eigencomponent and number of standard deviations
	SampleC<SArray1dC<Point2dC> > proj_pts = PerformComponentProjection(pca_mean, pca_sdev, pca_loadings, eigencomponent, devs, numcomponents);
	cout<<"Proj Pts "<<proj_pts<<endl;
	UIntT count = 0;
	
	ImageC<RealRGBValueC> img(IndexRange2dC(320,450,300,450));
	img.Fill(RealRGBValueC(255.0,255.0,255.0));
	for(SampleIterC<SArray1dC<Point2dC> > it(proj_pts); it; it++)
	{
		SArray1dC<Point2dC> pts = aff_mean + (*it);
		DisplayImage(pts,img,3);
	}
	if(!Save("@XA:Output Image",img)) cerr<<"Could not save the resulting output image"<<endl;
	*/
	//Check the PCA Projection and Affine Projection Classes
	SimilarityProjectionC sim(pca);
	
	SampleC<SArray1dC<Point2dC> > lip_feats;
	IStreamC lfile(lip);
	lfile>>lip_feats;
	//Now we have the lip features file - perform Affine projection
	SArray1dC<Point2dC> samplelip = lip_feats[0];
	cout<<"Projecting to affine space - "<<endl;
	Tuple2C<SArray1dC<Point2dC>, VectorC> aff = sim.ProjectToAffineSpace(samplelip);
	cout<<"Projected on to affine mean - "<<aff.Data1()<<"\n with resulting vector as "<<aff.Data2()<<endl;
	cout<<"Now projecting from affine space - "<<endl;
	SArray1dC<Point2dC> samplelipfromaffine = sim.ProjectFromAffineSpace(PointsToVector(aff.Data1()),aff.Data2());
	cout<<"The residual of the reverse projection is "<<SubtractSArrays(samplelipfromaffine,samplelip)<<endl;
	cout<<"About to perform PCA Projection of the residual data from the affine projection "<<endl;
	PCAProjectionC pcaproj(pca,numcomponents);
	cout<<"The input mean is "<<pcaproj.GetMean()<<endl;
	VectorC proj_vect = pcaproj.PerformPCAProjection(aff.Data1());
	cout<<"Projected vector in PCA Space is - "<<proj_vect<<endl;
	SArray1dC<Point2dC> inv_vect = pcaproj.PerformInverseProjection(proj_vect);
	cout<<"Projected vector from PCA Space is - "<<inv_vect<<endl;
	cout<<"Residual of the projection is "<<SubtractSArrays(inv_vect,aff.Data1())<<endl;
	cout<<"About to test ICA projection class"<<endl;
	ICAProjectionC ica(icadir);
	cout<<"Points to be input are : "<<PointsToVector(aff.Data1())<<endl;
	VectorC icavect = ica.PerformICAProjection(aff.Data1());
	cout<<"ICA Points are : "<<icavect<<endl;
	SArray1dC<Point2dC> origpts = ica.PerformInverseProjection(icavect);
	cout<<"Returned Points Are : "<<PointsToVector(origpts)<<endl;
	cout<<"Residual of projection is "<<SubtractSArrays(origpts,aff.Data1())<<endl;
	return 0;	
}

VectorC PerformPCAProjection(VectorC const &mu, VectorC const &lamda, MatrixC const &e)
{
	FuncMeanProjectionC pca(mu,e);
	return pca.Apply(lamda).Copy();
}

SampleC<SArray1dC<Point2dC> > PerformComponentProjection(VectorC const &mu, VectorC const &vals, MatrixC const &e, UIntT const &comp, IntT const &dev, UIntT const &ncp)
{
	cout<<"We are going to perform a projection along the "<<comp<<"th principal component with "<<dev<<" standard deviations and whilst retaining "<<ncp<<" dimensions"<<endl;
	VectorC lamda(ncp);
	lamda.Fill(0.0);
	SampleC<SArray1dC<Point2dC> > pts;
	for(IntT i = -(IntT)(dev); i <= (IntT)dev; i++)
	{
		lamda[comp] = (RealT)(vals[comp]*(RealT)i);
		//cout<<"lamda to be inserted is "<<lamda<<endl;
		MatrixC submat = (const_cast<MatrixC &>(e));
		VectorC res = PerformPCAProjection(mu,lamda,submat.SubMatrix(e.Rows(),ncp));
		pts.Append(VectorToSArray(res.Copy()).Copy());
	}
	return pts;	
}

SArray1dC<Point2dC> VectorToSArray(VectorC const &vec)
{
	// This is the case of the new R based representation where each vector is r1,c1,r1,c1
	UIntT size = vec.Size() / 2;
	SArray1dC<Point2dC> res(size);
	for(UIntT i = 0; i < size; i+=2)
	{
		Point2dC pt(vec[i],vec[i+1]);
		res[i] = pt.Copy();
	
	}
	return res;	
}
void DisplayImage(const SArray1dC<Point2dC> &pts, ImageC<RealRGBValueC> &img, const UIntT &size)
{
	//ImageC<RealRGBValueC> img(IndexRange2dC(0,700,0,700));
	//img.Fill(RealRGBValueC(0.0,0.0,0.0));
	//Generate Random Colour
	RealRGBValueC rgb = GenerateRandomColour();
	for(SArray1dIterC<Point2dC> it(pts); it; it++)
	{
		//if(img.Frame().Contains(Index2dC((*it).Row(),(*it).Col())))
			DrawCross(img,rgb,Index2dC((*it).Row(),(*it).Col()),size);
	}
	//~ if(!Save("@XA:Output Image",img)) cerr<<"Could not save the resulting output image"<<endl;
	//~ sleep(5);
}

RealRGBValueC GenerateRandomColour(void)
{
	//Generate 3 Random Numbers Between 1 and 255
	RealRGBValueC col((RealT)(myrandgen.Generate()*255),(RealT)(myrandgen.Generate()*255),(RealT)(myrandgen.Generate()*255));
	return col;
}

