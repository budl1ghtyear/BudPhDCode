
#ifndef UTILITYFUNCTIONS_hh
#define UTILITYFUNCTIONS_hh

//HEADER FILE TO STORE ALL THE UTILITY FUNCTIONS
#include "Ravl/Array2d.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/SArray1d.hh"
#include "Ravl/SArray1dIter.hh"
#include "Ravl/OS/Directory.hh"
#include "Ravl/OS/Filename.hh"
#include "Ravl/Matrix.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/Vector.hh"

using namespace RavlN;

template <typename T>
bool LoadFileFromStream(FilenameC const &file, T &dat)
//Load Some Data Type From An Input Stream Object
{
	IStreamC myfile(file);
	if(myfile.fail())
	{
		cerr<<"Could not find the filename - "<<file<<endl;
		return false;
	}
	else
	{
		myfile>>dat;
		myfile.Close();
		return true;
	}
}

template <typename T>
SArray1dC<T> ArrayToSArray_B(Array1dC<T> const &pts)
{
	SArray1dC<T> spts(pts.Size());
	for(UIntT i = 0; i < pts.Size(); i++)
	{
		spts[i] = pts[i].Copy();
	}
	return spts;
}

template <typename T>
SArray1dC<T> ArrayToSArray_S(Array1dC<T> const &pts)
{
	SArray1dC<T> spts(pts.Size());
	for(UIntT i = 0; i < pts.Size(); i++)
	{
		spts[i] = pts[i];
	}
	return spts;
}

template <typename T>
Array1dC<T> SArrayToArray_B(SArray1dC<T> const &pts)
{
	Array1dC<T> spts(pts.Size());
	for(UIntT i = 0; i < pts.Size(); i++)
	{
		spts[i] = pts[i].Copy();
	}
	return spts;
}

template <typename T>
SArray1dC<T> DListToSArray_B(DListC<T> const &pts)
{
	SArray1dC<T> sarray(pts.Size());
	UIntT i = 0;
	for(DLIterC<T> it(pts); it; it++)
		sarray[i++] = ((*it).Copy());
	return sarray;
}


template <typename T>
SArray1dC<T> DListToSArray_S(DListC<T> const &pts)
{
	SArray1dC<T> sarray(pts.Size());
	UIntT i = 0;
	for(DLIterC<T> it(pts); it; it++)
		sarray[i++] = (*it);
	return sarray;
}
template <typename T>
Array1dC<T> SArrayToArray_S(SArray1dC<T> const &pts)
{
	Array1dC<T> spts(pts.Size());
	for(UIntT i = 0; i < pts.Size(); i++)
	{
		spts[i] = pts[i];
	}
	return spts;
}

VectorC PointsToVector(SArray1dC<Point2dC> const &pts);

VectorC PointsToVector(Array1dC<Point2dC> const &pts);

SArray1dC<Point2dC> VectorToPoints(VectorC const &vect);

template <typename T>
SArray1dC<T> SubtractSArrays(SArray1dC<T> const &one, SArray1dC<T> const &two)
{
	if(one.Size()!= two.Size())
	{
		cerr<<"Two arrays are of unequal size"<<endl;
		exit(1);
	}
	SArray1dC<Point2dC> res(one.Size());
	for(UIntT i = 0; i < one.Size(); i++)
	{
		res[i] = one[i] - two[i];
	}
	return res;
}
template <typename T>
SArray1dC<T> AddSArrays(SArray1dC<T> const &one, SArray1dC<T> const &two)
{
	if(one.Size()!= two.Size())
	{
		cerr<<"Two arrays are of unequal size"<<endl;
		exit(1);
	}
	SArray1dC<Point2dC> res(one.Size());
	for(UIntT i = 0; i < one.Size(); i++)
	{
		res[i] = one[i] + two[i];
	}
	return res;
}

Point2dC ComputeCentroid(SArray1dC<Point2dC> const &pts);

#endif
