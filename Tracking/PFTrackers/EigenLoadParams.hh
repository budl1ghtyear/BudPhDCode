//EigenLoadParams.hh
//Author - Bud Goswami
//Purpose - to create a header file that contains the function definitions required to load the necessary eigen transform and AR estimation parameters
#ifndef EigenLoadParams_HH
#define EigenLoadParams_HH

#include 	"Ravl/IO.hh"
#include 	"Ravl/Matrix.hh"
#include 	"Ravl/OS/Filename.hh"

using namespace RavlN;

//Return a matrix by loading a textfile that does not specify size in the first line
//Input arguments - filename and the number of rows and cols you want
MatrixC GetMatrixNoSize(const FilenameC &f, const UIntT &r, const UIntT &c)
{
	IStreamC f_str(f);
	MatrixC out(r,c);
	for(UIntT i = 0; i < r; i++)
	{
		for(UIntT j = 0; j < c; j++)
		{
			f_str >> out[i][j];
		}
	}
	return out;
}
//Return a matrix by loading a textfile that specifies the size of the matrix in the first line
MatrixC GetMatrixWithSize(const FilenameC &f)
{
	IStreamC a1_str(f);
	IntT r = 0;
	IntT c = 0;
	a1_str>>r;
	a1_str>>c;
	//cout<<"R "<<r<<"\t C "<<c<<endl;
	MatrixC a1(r,c);
	for(IntT i = 0; i < r; i++)
	{
		for(IntT j = 0; j< c; j++)
		{
			a1_str>>a1[i][j];
		}
	}
	return a1;
}

#endif

