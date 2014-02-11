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
#include "MBspline.hh"
#include "ModellingUtilityFunctions.hh"
#include "Align.hh"
#include <cstdio>
#include "ASMRotationNormalisation.hh"
#include "ASMAffineNormalisation.hh"
#include "Ravl/PatternRec/SampleIter.hh"
#include "Ravl/Image/AAMAppearance.hh"
using namespace RavlN;
using namespace RavlImageN;

SampleC<AAMAppearanceC> LoadPointsFromFile(FilenameC const &file);

int main(int nargs,char **argv) 
{
	//CREATE THE REQUIRED INPUT OPTIONS
	OptionC opt(nargs,argv);
	//IO Options
	DirectoryC lip = opt.String("l","in.txt","Text File with Lip-Co-ordinates");
	StringC ext = opt.String("e","spl","Extension of the file containing the lip-coordinates in SampleC<SArray1dC<Point2dC> > format");
	RealT var = opt.Real("v",1.0,"Variation to be preserved in the parameter");
	IntT type = opt.Int("t",0,"Type of normalisation to be performed:0 - Perform 3 parameter GPA, 1 - Perform 4 parameter GPA, 2 - Perform 6 parameter GPA, 3 - Perform No GPA");
	opt.Check();
	
	//Load all the .aam files from an input directory
	DListC<StringC> aamfiles = lip.FiltList(StringC("*."+ext));
	//Now, for each file, load the aam data and append it to a global container
	ModellingUtilityFunctions mufc;
	SampleC<AAMAppearanceC> lip_feats;
	for(DLIterC<StringC> it(aamfiles); it; it++)
	{
		FilenameC file = lip + (*it);
		SampleC<AAMAppearanceC> pts = LoadPointsFromFile(file);
		lip_feats.Append(pts);
	}
	//Now we have a list of points in a SampleC<AppearanceC>
	SampleVectorC smpl; SArray1dC<Point2dC> mean;
	//Perform Generalised Procrustes Alignment Using Adapted Ravl AAM Classes
	switch(type)
	{
		case 0:
		{
			ASMNormalisationC asmdl(true);
			asmdl.Design(lip_feats, var,22);
			smpl = asmdl.GetNormalisedPoints();
			mean = asmdl.MeanPoints();
			break;
		}
		case 1:
		{	
			ASMScaleRotationNormalisationC asmdl(true);
			asmdl.Design(lip_feats, var,22);
			smpl = asmdl.GetNormalisedPoints();
			mean = asmdl.MeanPoints();
			break;
		}
		case 2:
		{
			ASMAffineNormalisationC asmdl(true);
			asmdl.Design(lip_feats, var,22);
			smpl = asmdl.GetNormalisedPoints();
			mean = asmdl.MeanPoints();
			break;
		}
		case 3:
		{
			SampleC<VectorC> res;
			UIntT dim = lip_feats.First().Points().Size() * 2;
			for(SampleIterC<AAMAppearanceC> it(lip_feats); it; it++)
			{
				VectorC vec(dim);
				UIntT count = 0;
				for(SArray1dIterC<Point2dC> sit((*it).Points()); sit; sit++)
				{
					vec[count] = (*sit).Row(); count++;
					vec[count] = (*sit).Col(); count++;
				}
				res.Append(vec.Copy());
			}
			smpl = res.Copy();
		}
	}
	//Output the resulting vector onto a screen
	for(SampleIterC<VectorC> it(smpl); it; it++)
	{
		UIntT size = it->Size();
		for(UIntT i = 0; i < size; i++)
		{
			cout<<(*it)[i]<<"\t";
		}
		cout<<endl;
	}
	//Save the mean points into a file
	FilenameC lip_file = lip+"MeanPoints.txt";
	OStreamC myfile(lip_file);
	myfile<<mean;
	return 0;	
}

SampleC<AAMAppearanceC> LoadPointsFromFile(FilenameC const &file)
{
	IStreamC myfile(file);
	if(myfile.fail())
	{
		cout<<"Could not open the file "<<file<<endl;
		exit(1);
	}	
	SampleC<SArray1dC<Point2dC> > arr;
	myfile>>arr;
	myfile.Close();
	SampleC<AAMAppearanceC> res;
	for(SampleIterC<SArray1dC<Point2dC> > it(arr); it; it++)
	{
		AAMAppearanceC aam((*it));
		res.Append(aam);
	}
	return res;
}

