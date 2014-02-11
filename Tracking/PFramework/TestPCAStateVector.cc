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
#include "PCAStateVector.hh"
#include "Ravl/RandomGauss.hh"
#include "ZOAVPropModel.hh"
#include "Particle.hh"
using namespace RavlN;
using namespace RavlImageN;
RandomGaussC BasePropagationModelC::rnd;
int main(int nargs,char **argv) 
{
	//CREATE THE REQUIRED INPUT OPTIONS
	OptionC opt(nargs,argv);
	//IO Options
	DirectoryC pca = opt.String("p","./PCA/","Directory containing the PCA data");
	UIntT numcomponents = opt.Int("n",10,"Number of principal components to keep");
	opt.Check();
	cout<<"Instantiate Sample Points"<<endl;
	SArray1dC<Point2dC> samplepts(11);
	samplepts[0] = Point2dC(379.957 , 336.994);
	samplepts[1] = Point2dC(377.137 , 343.896);
	samplepts[2] = Point2dC(364.833 , 355.647 );
	samplepts[3] = Point2dC(368.033 , 374.19 );
	samplepts[4] = Point2dC(363.711 , 394.17 );
	samplepts[5] = Point2dC(372.219 , 406.287 );
	samplepts[6] = Point2dC(381.719 , 419.781 );
	samplepts[7] = Point2dC(396.891 , 403.21);
	samplepts[8] = Point2dC(405.472 , 382.169 );
	samplepts[9] = Point2dC(399.442 , 352.066 );
	samplepts[10] = Point2dC(390.992 , 346.028 );
	cout<<"Input points are "<<samplepts<<endl;

	//Now we have the sample points instantiate the statevector
	PCAStateVectorC sv(pca,samplepts,numcomponents);
	//~ cout<<"Input Mean Is : "<<sv.GetPCA().GetMean()<<endl;
	//~ cout<<"Input Projection Is : "<<sv.GetPCA().GetEigenVectors()<<endl;
	//~ cout<<"Affine Mean Is : "<<sv.GetSimilarityProjection().GetMean()<<endl;
	cout<<"Resulting State Vector Is : "<<sv.GetStateVector()<<endl;
	cout<<"Resulting Back Projection Is: "<<sv.GetObservationVector()<<endl;
	cout<<"Resulting Residual Is : "<<SubtractSArrays(ArrayToSArray_B(sv.GetObservationVector()),samplepts)<<endl;
	PCAStateVectorC nsv(pca,sv.GetStateVector(),numcomponents);
	cout<<"Resulting Residual Is : "<<SubtractSArrays(ArrayToSArray_B(nsv.GetObservationVector()),samplepts)<<endl;
	//Now instantiate the particleC
	ParticleC pt(sv,1.0);
	ZOAVPropModelC pm(3.0);
	for(UIntT i = 0; i < 20; i++)
	{
		cout<<pm.Propagate(pt).GetState().GetStateVector()<<endl;
	}
	return 0;
}

