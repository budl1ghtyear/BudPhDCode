////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SampleData.cc
// Author: Bud Goswami (B.Goswami@surrey.ac.uk)
// Date: 25.02.08
//Defines the SampleData Class to store data that can have robust clustering performed on it
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Required Libraries
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SampleData.hh"

SampleData::SampleData()
{
	dims = 0;
	n = 0;
}

SampleData::~SampleData()
{}

SampleData::SampleData(const ImageC<RealRGBValueC> &inpimg,const RealRGBValueC &lbl)
{
	n = 0;
	for(Array2dIterC<RealRGBValueC> it(inpimg); it; it++)
	{
		if((*it)!=lbl)
		{
			RealRGBValueC rgb = *it;
			Index2dC curind = it.Index();
			VectorC vec(rgb.Red(), rgb.Green(), rgb.Blue(), curind.Row(), curind.Col());
			smpl.Append(vec);
			n++;
		}
	}
	dims = inpimg[inpimg.Frame().TopLeft()].Size();
}

SampleData::SampleData(const ImageC<VectorC> &inpimg)
{
	n = inpimg.Frame().Area();
	dims = inpimg[inpimg.Frame().TopLeft()].Size();
	for(Array2dIterC<VectorC> it(inpimg); it; it++)
	{
		VectorC pix = *it;
		Index2dC curind = it.Index();
		VectorC vec(pix[0],pix[1],pix[2],curind.Row(), curind.Col());
		smpl.Append(vec);
	}
}

SampleData::SampleData(const ImageC<TFVectorC<RealT,3> > &img)
{
	n = img.Frame().Area();
	dims = img[img.Frame().TopLeft()].Size();
	for(Array2dIterC<TFVectorC<RealT,3> > it(img); it ; it++)
	{
		VectorC pix(dims);
		for(IntT i = 0; i < dims; i++)
		{
			pix[i] = (*it)[i];
		}
		SArray1dC<RealT> indices(2);
		indices[0] = (RealT)(it.Index().Row()); indices[1] = (RealT)(it.Index().Col());
		pix.Append(indices);
		smpl.Append(pix);
	}
}
DListC<VectorC> SampleData::GetSmple()
{
	return smpl;
}

IntT SampleData::GetDimensions()
{
	return dims;
}
	
IntT SampleData::GetSize()
{
	return n;
}	
IntT SampleData::SetDimensions(const IntT &d)
{
	dims = d;
	return dims;
}

