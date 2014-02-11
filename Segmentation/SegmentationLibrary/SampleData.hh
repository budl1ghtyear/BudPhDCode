//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SampleData.hh
// Author: Bud Goswami (B.Goswami@surrey.ac.uk)
// Date: 25.02.08
//Declares the SampleData Class to store data that can have robust clustering performed on it
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SAMPLEDATA_H
#define SAMPLEDATA_H
//////////////////////////////////////////////////////////////////////////////
//! file =      "Ravl/LipCode/Segmentation/SegmentationLibrary/SampleData.hh"   
//! author =    "Bud Goswami"
//! lib =       SegLib
//! date =      "10.09.09"
//! userlevel = Normal
//! docentry =  "Ravl.API.LipCode.Segmentation.SegmentationLibrary"
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Required Libraries
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Ravl/Array2dIter.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/Index2d.hh"
#include "Ravl/TFVector.hh"
#include "Ravl/Vector.hh"
//Pixel Types
#include "Ravl/Image/RealRGBValue.hh"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace RavlN ;
using namespace RavlImageN;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//:Data structure required by the MCD methods
class SampleData
{
	public:
	//constructor functions
	SampleData();
	~SampleData();
	SampleData(const ImageC<RealRGBValueC> &inpimg, const RealRGBValueC &lbl);
	//: Constructor with a background labelled (with RealRGBValueC = lbl) ImageC<RealRGBValueC>
	SampleData(const ImageC<VectorC> &inpimg);
	//: Constructor with a vector image
	SampleData(const ImageC<TFVectorC<RealT,3> > &img);
	
	DListC<VectorC> GetSmple();
	IntT GetDimensions();
	IntT GetSize();
	IntT SetDimensions(const IntT &d);
	
	void operator=(const SampleData &other)
	{
		this->smpl = other.smpl.Copy();
		this->dims = other.dims;
		this->n = other.n;
	}
	
	private:
	DListC<VectorC> smpl; //: Store the data as a DList (this was used because the MCD method requires rapid sorting algorithm. DList implements MergeSort O(nLogn))
	IntT dims; //:store the number of dimensions
	IntT n; //:store size of sample
	
};

#endif
