//ClosedBSplineFitting Program
#include 	"Ravl/Array1d.hh"
#include 	"Ravl/Array1dIter.hh"
#include 	"Ravl/Circle2d.hh" 
#include 	"Ravl/DArray1d.hh"
#include 	"Ravl/DArray1dIter.hh"
#include 	"Ravl/Image/DrawCross.hh"
#include 	"Ravl/Image/Image.hh"
#include 	"Ravl/Image/RealRGBValue.hh"
#include 	"Ravl/IndexRange2d.hh"
#include 	"Ravl/IO.hh"
#include 	"Ravl/Option.hh"
#include 	"Ravl/Point2d.hh"
#include 	"Ravl/RandomGauss.hh"
#include 	"Ravl/RandomMersenneTwister.hh"

using namespace RavlN;
using namespace RavlImageN;
using namespace RavlConstN;
RandomGaussC myrandgen;
//Simple Functions to Perform Some Data Generation
Array1dC<Point2dC> GenerateCirclePoints(const UIntT &radius, const UIntT &npts); //Remember that the centre of the circle is at 50,50
void DisplayImage(const Array1dC<Point2dC> &pts, const UIntT &size=2);
RealRGBValueC GenerateRandomColour(void);

//Perform Knot Vector Generation 
Array1dC<RealT> PeriodicKnots(const UIntT &degree, const UIntT &nplusone);
Array1dC<RealT> Parameterisation(const Array1dC<Point2dC> &pts);
Array1dC<RealT> ComputeBasisFunction(const UIntT &c, const RealT &t, const UIntT &npts, const Array1dC<RealT> &x);
Array1dC<RealT> GetNonZeroVals(const IntT &i, const RealT &uval, const IntT &p, const Array1dC <RealT> &knots );
IntT GetKnotSpan(const IntT &n, const IntT &p, const RealT &uval, const Array1dC <RealT> &knots);
int main(int nargs,char **argv) 
{
	OptionC opt(nargs,argv);
	UIntT rad = opt.Int("r",100,"Circle Radius");
	UIntT npt = opt.Int("m",15,"Number of data Points to Generate");
	UIntT sze = opt.Int("s",5,"Size of display cross");
	UIntT deg =  opt.Int("d",3,"Degree of the B-Spline");
	UIntT npone = opt.Int("n",11,"Number of B-Spline control points (n+1)"); 
	opt.Check();
	Array1dC<Point2dC> circ_pts = GenerateCirclePoints(rad,npt);
	cout<<"Generated Points are - "<<circ_pts<<endl;
	DisplayImage(circ_pts,sze);
	Array1dC<RealT> kts = PeriodicKnots(deg,npone);
	cout<<"Generated Periodic Knots are - "<<kts<<endl;
	Array1dC<RealT> pars = Parameterisation(circ_pts);
	cout<<"Generated Chord Length Parameters are - "<<pars<<endl;
	/*for(Array1dIterC<RealT> it(pars); it; it++)
	{
		Array1dC<RealT> bas = ComputeBasisFunction((deg+1), (*it), npone, kts);
		//cout<<"Basis Function for parameter value "<<(*it)<<"\n"<<bas<<endl;
	}*/
	for(Array1dIterC<RealT> it(pars); it; it++)
	{
		Array1dC<RealT> bas = GetNonZeroVals(GetKnotSpan((IntT)npone,(IntT)deg+1,(*it),kts),(*it),deg+1,kts);
		cout<<"Basis Function for parameter value "<<(*it)<<"\n"<<bas<<endl;
	}	
	
	return 0;
}

Array1dC<Point2dC> GenerateCirclePoints(const UIntT &radius, const UIntT &npts)
{
	Array1dC<Point2dC> pts(npts); 
	Point2dC center(350,350);
	Circle2dC circle(center,radius);
	RealT incr = (2.0*pi)/(RealT)npts;
	for(UIntT i = 0; i < npts; i++)
	{
		RealT angle = (RealT)i * incr;
		Point2dC pt = circle.Value(angle);
		pts[i] = pt.Copy();
	}
	return pts;
}

void DisplayImage(const Array1dC<Point2dC> &pts, const UIntT &size)
{
	ImageC<RealRGBValueC> img(IndexRange2dC(0,700,0,700));
	img.Fill(RealRGBValueC(0.0,0.0,0.0));
	//Generate Random Colour
	RealRGBValueC rgb = GenerateRandomColour();
	for(Array1dIterC<Point2dC> it(pts); it; it++)
	{
		DrawCross(img,rgb,Index2dC((*it).Row(),(*it).Col()),size);
	}
	if(!Save("@XA:Output Image",img)) cerr<<"Could not save the resulting output image"<<endl;
}

RealRGBValueC GenerateRandomColour(void)
{
	//Generate 3 Random Numbers Between 1 and 255
	RealRGBValueC col((RealT)(myrandgen.Generate()*255),(RealT)(myrandgen.Generate()*255),(RealT)(myrandgen.Generate()*255));
	return col;
}

Array1dC<RealT> PeriodicKnots(const UIntT &degree, const UIntT &nplusone)
{
		Array1dC<RealT> knots(nplusone + degree); knots.Fill(0.0);
		for (UIntT i = 0 ; i < knots.Size(); i++)
		{
			knots[i] = (RealT)((RealT)(i) - (RealT)(degree)) / (RealT) (nplusone - degree); 
		}
		return knots;
}

Array1dC<RealT> Parameterisation(const Array1dC<Point2dC> &pts)
{
	Array1dC<RealT> pars(pts.Size()); pars.Fill(0.0); RealT chordlength = 0.0;
	for(UIntT i = 1; i < pars.Size(); i++)
	{
		chordlength += pts[i].EuclidDistance(pts[i-1]);
		pars[i] = chordlength;
	}
	for(UIntT i = 0; i < pars.Size(); i++)
	{
		pars[i] /= chordlength;
	}
	return pars;
}

Array1dC<RealT> ComputeBasisFunction(const UIntT &c, const RealT &t, const UIntT &npts, const Array1dC<RealT> &x)
//~ Algorithm from Intro to NURBS: 
//~ c = order of the B-Spline basis function
//~ d = first term of the basis function recursion solution
//~ e = second term of the basis function recursion solution
//~ npts = number of control polygon vertices
//~ n(,) = array containing the basis functions - n(1,1) contains the basis function associated with B1 etc
//~ nplusc = constant npts + c - maximum number of knot values
//~ t = parameter value
//~ x() = knot vector
{
	UIntT nplusc = npts + c;
	Array1dC<RealT> temp(IndexRangeC(1,nplusc-1)); temp.Fill(0.0);// allows for 20 polygon vertices
	//calculate the first order functions Ni,1
	for(UIntT i = 1; i <= nplusc - 1; i++)
	{
		if((t >= x[i]) && (t < x[i+1]))
			temp[i] = 1.0;
	}
	//calculate the higher order basis functions
	for(UIntT k = 2; k <= c; k++ )
	{
		RealT d = 0.0, e = 0.0;
		for(UIntT i = 1; i <= nplusc-k; i++)
		{
			if(temp[i] != 0.0)
			{
				if((x[i+k-1] - x[i]) != 0.0)
					d = ((RealT)(t - x[i])*temp[i])/(RealT)(x[i+k-1] - x[i]);
				else
					d = 0.0;
			}
			else
				d = 0.0;
			if(temp[i+1] != 0.0)
			{
				if((x[i+k] - x[i+1]) != 0.0)
					e = ((RealT)(x[i+k]-t)*(temp[i+1]))/((RealT)(x[i+k] - x[i+1]));
				else
					e = 0.0;
			}
			else
				e = 0.0;
			temp[i] = d+e;
		}
	}
	if (t == x[nplusc])
		temp[npts] = 1.0;
		
	Array1dC<RealT> n(IndexRangeC(1,npts));
	for(UIntT i = 1; i <= npts; i++)
	{
		n[i] = temp[i];
	}
	if(t == x[nplusc])
		n[npts] = 1.0;
		
	return n;
}

Array1dC<RealT> GetNonZeroVals(const IntT &i, const RealT &uval, const IntT &p, const Array1dC <RealT> &knots )
{
	IntT newi = i;
	Array1dC<RealT> vals(p+1);
	Array1dC<RealT> left(p+1);
	Array1dC<RealT> right(p+1);
	//fill the arrays with zeroes
	vals.Fill(0.0); left.Fill(0.0); right.Fill(0.0);
	vals[0] = 1.0;
	RealT saved ;
	RealT temp=0.0;
	for(IntT j = 1; j <= p; j++)
	{
		left[j] = uval - knots[newi+1-j];
		right[j] = knots[newi+j] - uval;
		saved = 0.0;
		for(IntT r = 0; r < j ; r++)
		{
			temp = vals[r]/(right[r+1] + left[j-r]);
			vals[r] = saved + right[r+1] * temp;
			saved = left[j-r] * temp;
			
		}
		vals[j] = saved;
	}
	cout<<"VALS - "<<vals<<endl;
	return vals;
}

IntT GetKnotSpan(const IntT &n, const IntT &p, const RealT &uval, const Array1dC <RealT> &knots)
{
	//first sort out the special case	
	if(uval == knots[n+1] )
		return n;
	IntT low = p;
	IntT high = n+1;
	//Do a binary search
	IntT mid = (low+high)/2;
	while( (uval < knots[mid] ) || (uval >= knots[mid+1] ) )
	{
		if (uval < knots[mid])
			high = mid;
		else 
			low = mid;
		mid = (low+high)/2;
	}
	cout<<"Knot Span for "<<uval<<" is "<<mid<<endl;
	return mid;
}
