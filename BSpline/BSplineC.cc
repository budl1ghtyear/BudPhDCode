#include "BSplineC.hh"

namespace RavlN {

//CONSTRUCTOR DEFINITIONS	
BSplineC::BSplineC(const UIntT &ord, const Array1dC<RealT> &kv,Array1dC<Point2dC> cp)
{
	order = ord;
	degree = ord - 1; //K-1
	numcp = cp.Size(); //N+1
	//open = true;
	knot_vec = kv.Copy(); //Knot Vector
	ctrl_vec = cp.Copy(); //Control Points
}

BSplineC::BSplineC(const UIntT &ord, const UIntT &ncp, KVG ktype, bool status)
{
	order = ord;
	degree = ord - 1; //K-1
	numcp = ncp;
	//numcp = ncp+degree;
	open = status;
	//if(!open)
	//{
	//	numcp += order;
	//}
	knot_vec = CalculateKnotVector(ord, ncp, ktype);	//To not have to play around with the rest of it
	//knot_vec = CalculateKnotVector(ord, ncp-ord, ktype);	//To not have to play around with the rest of it
}

///////////////////////////////////////////////////////////////////
//PROTECTED METHODS:
///////////////////////////////////////////////////////////////////

//PARAMETER SELECTION METHODS
//ChordLength() - as specified by the CS MTU notes
Array1dC<RealT> BSplineC::ChordLengthParSel(const Array1dC<Point2dC> &points)
{
	Array1dC<RealT> out(points.Size());
	Array1dC<RealT> cursum(points.Size());
	out.Fill(0.0); cursum.Fill(0.0);
	RealT totlen = 0.0; //totlen - contains the total chord length
	for(UIntT it = 1; it < out.Size(); it++)
	{
		Point2dC now = points[it].Copy();
		Point2dC nxt = points[it-1].Copy();
		RealT sum = Sqrt(Sqr(now.Row() - nxt.Row()) + Sqr(now.Col() - nxt.Col())); //calculate square euclid distance
		out[it] = sum;
		totlen += sum;
		cursum[it] = totlen;
	}
	cursum /= totlen;
	return cursum;
}	
//CentripetalParSel()
Array1dC<RealT> BSplineC::CentripetalParSel(const Array1dC<Point2dC> &points)
{
	Array1dC<RealT> out(points.Size());
	out.Fill(0.0);
	RealT a = 0.5; //a is the power value to which the distances must be raised
	RealT totlen = 0.0; //totlen - contains the total chord length
	for(UIntT it = 1; it < out.Size(); it++)
	{
		Point2dC now = points[it].Copy();
		Point2dC nxt = points[it-1].Copy();
		RealT sum = Pow(Sqrt(Sqr(now.Row() - nxt.Row()) + Sqr(now.Col() - nxt.Col())),a); //calculate square euclid distance
		totlen += sum;
		out[it] = totlen;
	}
	out /= totlen;
	return out;	
}
//Uniform Method()
Array1dC<RealT> BSplineC::UniformParSel(const Array1dC<Point2dC> &points)
{
	Array1dC<RealT> out(points.Size());
	out.Fill(0.0);
	RealT n = static_cast<RealT>(points.Size()-1); //to store the value of n (tot num of data points is n+1)
	for(IndexC i = 1; i < points.Size(); i++)
	{
		//if(i == points.Size()-1)
		//	out[i] = 1.0;
		//else
		{
			RealT val = (static_cast<RealT>(i))/n;
			out[i] = val;
		}
	}
	return out;	
}

//KNOT VECTOR GENERATION
//UniformPeriodic() e.g. [0123]
Array1dC<RealT> BSplineC::UniformPeriodicKVG(const UIntT &ord, const UIntT &nc)
{
	//For algorithm, see page 47 of Introduction to NURBS
	//Remember - number of knots = n+1+k (n+1 - number of control points, k = order)
	Array1dC<RealT> knots(IndexRangeC(1,(nc+ord))); //Remember that 1<= i <= n+1+k
	knots.Fill(0.0);
	for(UIntT i = 2; i <= (nc+ord); i++)
	{
		knots[i] = i - 1;
	}
	return knots;
}
//Uniform Open() e.g. [00012333]
Array1dC<RealT> BSplineC::UniformOpenKVG(const UIntT &ord, const UIntT &nc)
{
	//For algorithm, see page 47 of Introduction to NURBS
	//Remember - number of knots = n+1+k (n+1 - number of control points, k = order)
	Array1dC<RealT> knots(IndexRangeC(1,(nc+ord))); //Remember that 1<= i <= n+1+k
	knots.Fill(0.0);
	//Fill array
	for(UIntT i = 1; i <= (nc+ord); i++)
	{
		//For first few knots, default value is zero anyway, so no need for explicit implementation
		if(((ord+1)<= i)&&(i <= nc)) //x_i = i - k, if (k+1) <= i <= (n+1)
			knots[i] = i - ord;
		else if(((nc+1)<= i)&&(nc+ord)) // x = n-k+2, if (n+2)<= i <= (n+k+1)
			knots[i] = (nc - 1) - ord + 2;
	
	}
	return knots;
}
//Non-Uniform Periodic() 
Array1dC<RealT> BSplineC::NonUniformPeriodicKVG(const UIntT &ord, const UIntT &nc, const Array1dC<Point2dC> &cp)
{
	cout<<"METHOD NOT IMPLEMENTED YET - RETURNING UNIFORM PERIODIC KV"<<endl;
	Array1dC<RealT> out = UniformOpenKVG(ord, nc);
	return out;
}
//Non-Uniform Open() 
Array1dC<RealT> BSplineC::NonUniformOpenKVG(const UIntT &ord, const UIntT &nc, const Array1dC<Point2dC> &cp)
{
	Array1dC<RealT> out(cp.Size());
	Array1dC<RealT> cursum(cp.Size());
	out.Fill(0.0); cursum.Fill(0.0);
	RealT totlen = 0.0; //totlen - contains the total chord length
	for(UIntT it = 1; it < out.Size(); it++)
	{
		Point2dC now =cp[it].Copy();
		Point2dC nxt = cp[it-1].Copy();
		RealT sum = Sqrt(Sqr(now.Row() - nxt.Row()) + Sqr(now.Col() - nxt.Col())); //calculate square euclid distance
		out[it] = sum;
		totlen += sum; 
		cursum[it] = totlen;
	}
	//totlen now contains total chordlength
	//out now contains the chordlengths of the B-Spline control points
	//For algorithm, see page 66 of Introduction to NURBS
	//Remember - number of knots = n+1+k (n+1 - number of control points, k = order)
	Array1dC<RealT> knots(IndexRangeC(1,(nc+ord))); //Remember that 1<= i <= n+1+k
	knots.Fill(0.0);
	//Fill array
	for(UIntT i = 1; i <= (nc+ord); i++)
	{
		//For first few knots, default value is zero anyway, so no need for explicit implementation
		if((1 <= i)&&(i <= (nc - ord)))
		{
			knots[i+ord] = ((RealT(((RealT)(i)/(RealT)(nc - 1 - ord + 2))*(out[i+1])) + (cursum[i]))/(totlen)) * (RealT)((nc-1)-ord + 2);
		}
		else if(((nc+1)<= i)&&(nc+ord))
		{
		 knots[i] = (nc - 1) - ord + 2;
		}
	}
	return knots;
}
//Closed Uniform and Periodic Knot Vector
Array1dC<RealT> BSplineC::ClosedUniformPeriodicKVG(const UIntT &ord, const UIntT &nc)
{
	//First get uniform perdiodic KV
	Array1dC<RealT> kv = this->UniformPeriodicKVG(ord, nc);
	RealT kvSize = (RealT)kv.Size();
	//Now wrap the knots so that the first k and the last k knots are the same
	Array1dC<RealT> addedknotTop(ord/2);
	Array1dC<RealT> addedknotBottom(ord/2);
	addedknotTop.Fill((RealT)0.0);
	addedknotBottom.Fill((RealT)0.0);
	for(UIntT i = 0; i < addedknotTop.Size(); i++)
	{	addedknotTop[i] = -(RealT)(addedknotTop.Size()-i);
		addedknotBottom[i] = (RealT)(kv.Size()) + i;
	}
	//kv = kv.Append(addedknot);
	kv=addedknotTop.Append(kv).Copy();
	kv = kv.Append(addedknotBottom).Copy();
	kv /= (kvSize-1);
	cout<<"Appended Knot Vector Top = "<<addedknotTop<<endl;
	cout<<"Appended Knot Vector Bot = "<<addedknotBottom<<endl;
	cout<<"Closed Knot Vector = "<<kv<<endl;
	return kv;
}



//Get basis function values
Array1dC<RealT> BSplineC::ComputeBasis(const UIntT &ord, const UIntT &nc, const RealT &par, const Array1dC<RealT> &kv) const
{
	UIntT nplusc = nc + ord;
	RealT d = 0.0;
	RealT e = 0.0;
	Array1dC<RealT> temp(IndexRangeC(1,nplusc-1));
	//Compute 1st Order Values
	for(UIntT i = 1; i < temp.Size(); i++)
	{
		if((par >= kv[i])&&(par < kv[i+1]))
			temp[i] = 1;
		else 
			temp[i] = 0;
	}
	/* calculate the higher order basis functions */
	for(UIntT k = 2; k <= ord; k++)
	{
		for(UIntT i = 1; i <= nplusc - k; i++)
		{
			//Compute First Term of Recursion
			if (temp[i] != 0)    /* if the lower order basis function is zero skip the calculation */
			{	
				RealT den = kv[i+k-1] - kv[i];
				if(den != 0.0)
					d = RealT((RealT)((par - kv[i])*temp[i]) / den);
				else 
					d = 0.0;
			}
			else
				d = 0.0;
			
			//Compute Second Term of Recursion
			if (temp[i+1] != 0)     /* if the lower order basis function is zero skip the calculation */
			{
				RealT den = kv[i+k] - kv[i+1];
				if(den != 0)
					e = (RealT(kv[i+k] - par)*temp[i+1]) / den;
				else
					e = 0.0;
			}
			else
				e = 0.0;
			
			temp[i] = d+e;
		}
	}
	if (par ==  kv[nplusc]){		/*    pick up last point	*/
 		temp[nc] = 1;
	}
	Array1dC<RealT> res(nc);
	for(UIntT i = 0; i < nc; i++)
	{
		res[i] = temp[i+1];
	}
	return res;
}

Array1dC<RealT> BSplineC::ComputeFirstDer(const UIntT &ord, const UIntT &nc, const RealT &par, const Array1dC<RealT> &kv) const
{
	UIntT nplusc = nc + ord;
	RealT d = 0.0;
	RealT e = 0.0;
	RealT f1 = 0.0, f2 = 0.0, f3 = 0.0, f4 = 0.0;
	Array1dC<RealT> temp(IndexRangeC(1,nplusc-1)); temp.Fill(0.0);
	Array1dC<RealT> temp1(IndexRangeC(1,nplusc-1)); temp1.Fill(0.0);
	//Compute 1st Order Values
	for(UIntT i = 1; i < temp.Size(); i++)
	{
		if((par >= kv[i])&&(par < kv[i+1]))
			temp[i] = 1;
		else 
			temp[i] = 0;
	}
	if (par ==  kv[nplusc])
	{		/*    pick up last point	*/
 		temp[nc] = 1;
	}
	/* calculate the higher order basis functions */
	for(UIntT k = 2; k <= ord; k++)
	{
		for(UIntT i = 1; i <= nplusc - k; i++)
		{
			//Compute First Term of Recursion
			if (temp[i] != 0.0)    /* if the lower order basis function is zero skip the calculation */
			{	
				RealT den = kv[i+k-1] - kv[i];
				if(den != 0.0)
					d = RealT((RealT)((par - kv[i])*temp[i]) / den);
				else 
					d = 0.0;
			}
			else
				d = 0.0;
			
			//Compute Second Term of Recursion
			if (temp[i+1] != 0.0)     /* if the lower order basis function is zero skip the calculation */
			{
				RealT den = kv[i+k] - kv[i+1];
				if(den != 0.0)
					e = (RealT(kv[i+k] - par)*temp[i+1]) / den;
				else
					e = 0.0;
			}
			else
				e = 0.0;
			
			//Insert First Derivative Stuff in Here
			//First term
			if (temp[i] != 0.0)       /* if the lower order basis function is zero skip the calculation */
			{
				RealT den = kv[i+k-1] - kv[i];
				if(den != 0.0)
	            f1 = temp[i]/den;
	         else
	         	f1 = 0.0;
	       }
	       else
	            f1 = 0.0;
	      //Second Term with only basis term
        	if (temp[i+1] != 0.0)     /* if the lower order basis function is zero skip the calculation */
        	{
        		RealT den = kv[i+k] - kv[i+1];
        		if(den != 0.0)
        			f2 = -temp[i+1]/den;
        		else
        			f2 = 0.0;
        	}
	     	else
	      	f2 = 0;
			
			//Third term (first term with a derivative calculation)
			if (temp1[i] != 0.0)      /* if the lower order basis function is zero skip the calculation */
			{
				RealT den = kv[i+k-1] - kv[i];
				if(den != 0.0)
					f3 = ((par-kv[i])*temp1[i]) / den;
				else
					f3 = 0.0;
			}
	      else
	         f3 = 0;			
	      //Fourth term 
			if (temp1[i+1] != 0)    /* if the lower order basis function is zero skip the calculation */
	      {
	      	RealT den = kv[i+k] - kv[i+1];
	      	if(den != 0.0)
	      		f4 = ((kv[i+k]-par)*temp1[i+1])/den;
	      	else
	      		f4 = 0.0;
	      }
	      else
	        	f4 = 0;			
			
			temp[i] = d+e;
			temp1[i] = f1 + f2 + f3 + f4;
		}
	}
	// Please note that the above code uses the dbasis.c implementation (with the same indexing system)
	//t = par, x[] = kv[]
	Array1dC<RealT> res(nc);
	for(UIntT i = 0; i < nc; i++)
	{
		res[i] = temp1[i+1];
	}
	return res;
}


MatrixC BSplineC::GetBasisMatrix(const UIntT &ord, const UIntT &nc, const Array1dC<RealT> &pars, const Array1dC<RealT> &kv)
{
	//Normalise knot vector
	//cout << "Knot vector for basis matrix computation:" << kv << endl;
	Array1dC<RealT> kvalt(kv.Copy());
	RealT kmax = kvalt[kvalt.Size() - 1];
	for(Array1dIterC<RealT> it(kvalt); it; it++)
	{
		(*it) /= kmax;
	}
	MatrixC bas(pars.Size(),nc);
	
	UIntT rowind = 0;
	for(Array1dIterC<RealT> it(pars); it; it++)
	{
		/*if(it.Index() == pars.Size()-1)
		{
			Array1dC<RealT> bas_val = this->ComputeBasis(ord, nc, (*it),kvalt);
			bas_val.Fill(0.0);
			bas_val[bas_val.Size()-1] = 1.0; 
			UIntT colind = 0;
			for(Array1dIterC<RealT> it2(bas_val); it2; it2++)
			{
				bas[rowind][colind] = (*it2);
				colind++;
			}			
		}
		else*/
		{
			Array1dC<RealT> bas_val = this->ComputeBasis(ord, nc, (*it),kvalt);
			//Now given the basis function values, insert it into the appropriate row in the matrix
			UIntT colind = 0;
			for(Array1dIterC<RealT> it2(bas_val); it2; it2++)
			{
				bas[rowind][colind] = (*it2);
				colind++;
			}
		}
		rowind++;
	}	
	return bas;
}

//Normalise Knot Vector (required for most internal methods)
Array1dC<RealT> BSplineC::NormaliseKnotVector(const Array1dC<RealT> &kv)
{
	Array1dC<RealT> kvalt(kv.Copy());
	RealT kmax = kvalt[kvalt.Size() - 1];
	for(Array1dIterC<RealT> it(kvalt); it; it++)
	{
		(*it) /= kmax;
	}		
	return kvalt;
}

MatrixC BSplineC::GetMatrixFromArray(const Array1dC<Point2dC> &points)
{
	MatrixC dat_mat(points.Size(), 2);
	UIntT rowind = 0;
	for(Array1dIterC<Point2dC> it(points); it; it++)
	{
		dat_mat[rowind][0] = (*it).Row();
		dat_mat[rowind][1] = (*it).Col();
		rowind++;
	}
	return dat_mat;
}
///////////////////////////////////////////////////////////////////
//PUBLIC METHODS:
///////////////////////////////////////////////////////////////////

//GetParameters()
Array1dC<RealT> BSplineC::GetParameters(const Array1dC<Point2dC> &points, ParSel ptype)
{
	Array1dC<RealT> pars;
	switch(ptype)
	{
		case UNIFORM:
		{
			pars = UniformParSel(points).Copy();
			break;
		}
		case CHORDLENGTH:
		{
			pars = ChordLengthParSel(points).Copy();
			break;
		}
		case CENTRIPETAL:
		{
			pars = CentripetalParSel(points).Copy();
			break;
		}
	};
	return pars;
}
//Obtain Knot Vector
Array1dC<RealT> BSplineC::CalculateKnotVector(const UIntT &ord, const UIntT &nc, KVG ktype)
{
	Array1dC<RealT> knots;
	switch(ktype)
	{
		case UPERIODIC: 
		{
			knots = UniformPeriodicKVG(ord, nc).Copy();
			break;
		}
		case UOPEN: 
		{
			knots = UniformOpenKVG(ord, nc).Copy();
			break;
		}
		case NPERIODIC: 
		{
			cout<<"Cannot Use Non-uniform Periodic KV Gen without specifying control points"<<endl;
			break;
		}
		case NOPEN: 
		{
			cout<<"Cannot Use Non-Uniform Open KV Gen without specifying control points"<<endl;
			break;
		}	
		case UPCLOSED:
		{
			knots = ClosedUniformPeriodicKVG(ord, nc).Copy();
			break;
		}
	};
	this->SetKnotVector(knots);
	return knots;
}
//Different Implementation of the CalculateKnotVector using different input arguments
Array1dC<RealT> BSplineC::CalculateKnotVector(const UIntT &ord, const UIntT &nc, const Array1dC<Point2dC> &cp, KVG ktype)
{
	Array1dC<RealT> knots;
	switch(ktype)
	{
		case NPERIODIC: 
		{
			knots = this->NonUniformPeriodicKVG(ord, nc, cp).Copy();
			break;
		}
		case NOPEN: 
		{
			knots = this->NonUniformOpenKVG(ord, nc, cp).Copy();
			break;
		}	
		case UPERIODIC: 
		{
			knots = this->UniformPeriodicKVG(ord, nc).Copy();
			break;
		}
		case UOPEN: 
		{
			knots = this->UniformOpenKVG(ord, nc).Copy();
			break;
		}
		case UPCLOSED:
		{
			knots = ClosedUniformPeriodicKVG(ord, nc).Copy();
			break;
		}
	};
	this->SetKnotVector(knots);
	return knots;
}
//CalculateControlPoints() - Performs Curve Approximation Given Parameters
Array1dC<Point2dC> BSplineC::CalculateControlPoints(const Array1dC<Point2dC> &points, ParSel ptype, KVG ktype, const UIntT &ord, const UIntT &ncp)
{
	order = ord; numcp = ncp;
	//Given the data point, perform parameter selection depending on ParSel
	//Array1dC<RealT> pars;
	//if (!open)
	//{
	//	//numcp += order;
	//	Array1dC<Point2dC> pointsWrapped(points.Size()+order);
	//	for (Array1dIter2C<Point2dC,Point2dC> it(points,pointsWrapped,false);it;it++)
	//	{
	//		it.Data2() = it.Data1().Copy();
	//	}
	//	for (IndexC i=points.Size();i<pointsWrapped.Size();i++) pointsWrapped[i]=points[i-points.Size()].Copy();
	//	pars = this->GetParameters(pointsWrapped, ptype);
	//}
	//else
	//	pars = this->GetParameters(points, ptype);
	Array1dC<RealT> pars = this->GetParameters(points, ptype);

	cout<<"Parameters are - "<<pars<<endl;
	Array1dC<RealT> knots = this->CalculateKnotVector(ord, ncp, points, ktype);
	//Array1dC<RealT> knots = this->CalculateKnotVector(ord, ncp-ord, points, ktype);
	knot_vec = knots.Copy();
	cout<<"Knots are - "<<knot_vec<<endl;

	MatrixC bas_mat = this->GetBasisMatrix(ord, ncp, pars, knots);
	//MatrixC bas_mat = this->GetBasisMatrix(ord, ncp-ord, pars, knots);
	cout<<"Basis Matrix "<<bas_mat<<endl;
	//Now depending on whether or not least squares estimation needs to be performed, sort out solution
	//Get D matrix
	MatrixC dat_mat = this->GetMatrixFromArray(points);
	MatrixC sol;
	if(bas_mat.Rows() == bas_mat.Cols()) //if the system is square
		sol = Solve(bas_mat, dat_mat);
	else
	 	sol = Solve(bas_mat.ATA(),MatrixC(bas_mat.T()*dat_mat));
	//At this point, we should have a solution matrix with ncp rows and 2 cols, so put them back into an array
	Array1dC<Point2dC> out(ncp);
	UIntT rowind = 0;
	for(Array1dIterC<Point2dC> it(out); it; it++)
	{
		(*it) = Point2dC(sol[rowind][0], sol[rowind][1]);
		rowind++;
	}	
	ctrl_vec = out.Copy();
	return out;
} 

Array1dC<Point2dC> BSplineC::CalculateControlPoints(const Array1dC<Point2dC> &points, KVG ktype, ParSel ptype)
{
	Array1dC<Point2dC> cp = this->CalculateControlPoints(points, ptype, ktype, order, numcp);
	return cp;
}

Array1dC<Point2dC> BSplineC::CalculateControlPoints(const Array1dC<Point2dC> &points, const Array1dC<RealT> &kv, ParSel ptype)
{
	//Array1dC<RealT> pars;
	//if (!open)
	//{
	//	Array1dC<Point2dC> pointsWrapped(points.Size()+order);
	//	for (Array1dIter2C<Point2dC,Point2dC> it(points,pointsWrapped,false);it;it++)
	//	{
	//		it.Data2() = it.Data1().Copy();
	//	}
	//	for (IndexC i=points.Size();i<pointsWrapped.Size();i++) pointsWrapped[i]=points[i-points.Size()].Copy();
	//	pars = this->GetParameters(pointsWrapped, ptype);
	//}
	//else
	//	pars = this->GetParameters(points, ptype);
	Array1dC<RealT> pars = this->GetParameters(points, ptype);

	knot_vec = kv.Copy();	
	//cout<<"Knots are - "<<knot_vec<<endl;
	MatrixC bas_mat = this->GetBasisMatrix(order, numcp, pars, knot_vec);
	//MatrixC bas_mat = this->GetBasisMatrix(ord, ncp, pars, knots);
	//cout<<"Basis Matrix "<<bas_mat<<endl;
	//Now depending on whether or not least squares estimation needs to be performed, sort out solution
	//Get D matrix
	MatrixC dat_mat = this->GetMatrixFromArray(points);
	MatrixC sol;
	if(bas_mat.Rows() == bas_mat.Cols()) //if the system is square
		sol = Solve(bas_mat, dat_mat);
	else
	 	sol = Solve(bas_mat.ATA(),MatrixC(bas_mat.T()*dat_mat));
	//At this point, we should have a solution matrix with ncp rows and 2 cols, so put them back into an array
	Array1dC<Point2dC> out(numcp);
	UIntT rowind = 0;
	for(Array1dIterC<Point2dC> it(out); it; it++)
	{
		(*it) = Point2dC(sol[rowind][0], sol[rowind][1]);
		rowind++;
	}	
	ctrl_vec = out.Copy();
	return out;	
}

//Render the curve with the current parameters and knot_vector configuration
Array1dC<Point2dC> BSplineC::RenderCurve(const UIntT &res)
{
	Array1dC<RealT> kvalt = NormaliseKnotVector(knot_vec.Copy());
	//Create a new KV from Game Phy
	//Array1dC<RealT> kvalt = this->NormaliseKnotVector(UniformOpenKVG(order, numcp+1));
	//Array1dC<RealT> kvalt = knot_vec.Copy();
	//cout<<"Normalised Knot Vector "<<kvalt<<endl;
	Array1dC<RealT> pars(res);
	RealT inc = 1.0/(RealT)(res-1);
	RealT val = 0.0;
	for(Array1dIterC<RealT> it(pars); it; it++)
	{
		(*it) = val;
		val += inc;
	}
	//cout<<"Parameter values "<< pars << endl;
	//Now pars contains the parameter values to be drawn
	MatrixC cp_mat(ctrl_vec.Size(), 2);
	//MatrixC cp_mat(ctrl_vec.Size()+1, 2);
	UIntT rowind = 0;
	for(Array1dIterC<Point2dC> it(ctrl_vec); it; it++)
	{
		cp_mat[rowind][0] = (*it).Row();
		cp_mat[rowind][1] = (*it).Col();
		rowind++;
	}	
	//cp_mat[cp_mat.Rows()-1][0] = cp_mat[0][0]; //GPHY - SUGGESTS ADDING AN EXTRA CONTROL POINT
	//cp_mat[cp_mat.Rows()-1][1] = cp_mat[0][1]; //GPHY - SUGGESTS ADDING AN EXTRA CONTROL POINT
	//cout<<"Control Point Matrix = "<<cp_mat<<endl;
	Array1dC<Point2dC> out(res);
	rowind = 0;
	for(Array1dIterC<RealT> it(pars); it; it++)
	{
		Array1dC<RealT> bas = this->ComputeBasis(order, numcp, (*it),kvalt);
		MatrixC bas_mat(1,bas.Size());
		for(UIntT i = 0; i < bas.Size(); i++)
			bas_mat[0][i] = bas[i];		
		MatrixC data = bas_mat*cp_mat;
		out[rowind] = Point2dC(data[0][0], data[0][1]);
		rowind++;
	}	
	return out;
}
//Render a closed B-Spline curve using the wrapping control points technique
Array1dC<Point2dC> BSplineC::RenderClosedCurve(const UIntT &res)
{
	Array1dC<RealT> kvalt = knot_vec.Copy();
	//Normalise Knot Vector
	RealT kmax = kvalt[kvalt.Size() - 1];
	for(Array1dIterC<RealT> it(kvalt); it; it++)
	{
		(*it) /= kmax;
	}	
	//Wrap Knots
	//Array1dC<RealT> kv_alt(order + 1); //to store the new k+1 knots
	Array1dC<RealT> kv_alt(order); //to store the new k+1 knots
	kv_alt.Fill(0.0);
	//for(UIntT i = 0; i < kv_alt.Size(); i++)
	for(UIntT i = 0; i < order; i++)
	{
		kv_alt[i] = kvalt[i];
	}
	//cout<<"Addendum to Knot Vector = "<<kv_alt<<endl;
	Array1dC<RealT> closedkv = kvalt.Append(kv_alt);
	//cout<<"Closed Knot Vector = "<<closedkv<<endl;
	//Now we have a closed knot vector
	Array1dC<RealT> pars(res);
	RealT inc = 1.0/(RealT)(res);
	RealT val = 0.0;
	for(Array1dIterC<RealT> it(pars); it; it++)
	{
		(*it) = val;
		val += inc;
	}
	//Now pars contains the parameter values to be drawn (also remember that P_n+1 = P_0)
	MatrixC cp_mat(ctrl_vec.Size()+1, 2);
	UIntT rowind = 0;
	for(Array1dIterC<Point2dC> it(ctrl_vec); it; it++)
	{
		cp_mat[rowind][0] = (*it).Row();
		cp_mat[rowind][1] = (*it).Col();
		rowind++;
	}	
	cp_mat[rowind][0] = ctrl_vec[0].Row();
	cp_mat[rowind][1] = ctrl_vec[0].Col();
	//Control Point Matrix Now Obtained
	Array1dC<Point2dC> out(res);
	rowind = 0;
	for(Array1dIterC<RealT> it(pars); it; it++)
	{
		//Array1dC<RealT> bas = this->ComputeBasis(order, numcp, (*it),closedkv);
		Array1dC<RealT> bas = this->ComputeBasis(order, numcp+1, (*it),closedkv);
		MatrixC bas_mat(1,bas.Size());
		for(UIntT i = 0; i < bas.Size(); i++)
			bas_mat[0][i] = bas[i];		
		MatrixC data = bas_mat*cp_mat;
		out[rowind] = Point2dC(data[0][0], data[0][1]);
		rowind++;
	}	
	return out;
}
//GetNormalLine for a given parameter value
LinePP2dC BSplineC::GetNormalLine(const RealT &par) const
{
	//Array1dC<RealT> kvalt = NormaliseKnotVector(knot_vec.Copy());
	Array1dC<RealT> kvalt = knot_vec.Copy();	//ALREADY NORMALISED - no need for the above
	Array1dC<RealT> derbas = ComputeFirstDer(order, numcp,par, kvalt);
	MatrixC db_mat(1,derbas.Size());
	for(UIntT i = 0; i < derbas.Size(); i++)
	{
		db_mat[0][i] = derbas[i];
	}
	//Get CP Matrix
	MatrixC cp_mat(ctrl_vec.Size(), 2);
	UIntT rowind = 0;
	for(Array1dIterC<Point2dC> it(ctrl_vec); it; it++)
	{
		cp_mat[rowind][0] = (*it).Row();
		cp_mat[rowind][1] = (*it).Col();
		rowind++;
	}	
	//cout<<"DB Mat = "<<db_mat<<endl;
	MatrixC point = db_mat * cp_mat;
	//cout<<"Point is "<<point<<endl;
	//For Tangent Point (x',y') - normal point = (-y',x')
	//Vector2dC normal(-point[0][1],point[0][0]);	//Using Tangent Just to Check
	Vector2dC normal(-point[0][1],point[0][0]);
	normal = normal.MakeUnit();
	Point2dC curvept = this->GetPointOnCurve(par);
	//cout<<"For curve Point = "<<curvept<<" the normal point is = "<<normal<<endl;
	//cout<<"manual angle calc = "<<ATan2(normal[1],normal[0])<<endl;
	LinePP2dC out(curvept, normal);
	return out;
}
//Given a value for the parameter, obtain the point on the curret B-Spline Curve
Point2dC BSplineC::GetPointOnCurve(const RealT &par) const
{
	//Array1dC<RealT> kvalt = NormaliseKnotVector(knot_vec.Copy());
	Array1dC<RealT> kvalt = knot_vec.Copy();	//NO NEED TO NORMALISE - already done
	Array1dC<RealT> bas = this->ComputeBasis(order, numcp, par, kvalt);
	MatrixC bas_mat(1,bas.Size());
	for(UIntT i = 0; i < bas.Size(); i++)
			bas_mat[0][i] = bas[i];		
	//Get CP Matrix
	MatrixC cp_mat(ctrl_vec.Size(), 2);
	UIntT rowind = 0;
	for(Array1dIterC<Point2dC> it(ctrl_vec); it; it++)
	{
		cp_mat[rowind][0] = (*it).Row();
		cp_mat[rowind][1] = (*it).Col();
		rowind++;
	}	
	MatrixC point = bas_mat * cp_mat;
	Point2dC pt(point[0][0],point[0][1]);
	return pt;
}
}
