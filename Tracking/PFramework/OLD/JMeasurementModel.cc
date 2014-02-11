#include "JMeasurementModel.hh"
JMeasurementModel::JMeasurementModel(const MeanCovarianceC &lp, const MeanCovarianceC &sk, const ImageC<RealRGBValueC> &img, const IntT &wd)
{
	lip = lp.Copy();
	skin = sk.Copy();
	origimg = img.Copy();
	window = wd;
	cout<<"MModel Constr Called"<<endl;
}
ImageC<VectorC> JMeasurementModel::GetNormRGVals(const ImageC<RealRGBValueC> &im)
{
	ImageC<VectorC> res(im.Frame());
	VectorC zer(0.0,0.0);
	res.Fill(zer);
	for(Array2dIter2C<RealRGBValueC, VectorC> it(im, res); it; it++)
	{
		RealRGBValueC pix = it.Data1().Copy();
		RealT intensity = pix.Y()*3.0;
		VectorC vec( pix.Red() / intensity,pix.Green() / intensity);
		it.Data2() = vec.Copy();
	}
	return res;  	
}

RealT JMeasurementModel::Measure(ParticleC &pt)
{
	//Array1dC<Point2dC> ind = VectortoPointArray(pt.GetState().GetX());
	Array1dC<Point2dC> ind = pt.GetState().GetObservationVector();
	RealT weight = GetJWeight(ind);
	return weight;
}

Array1dC<Index2dC> JMeasurementModel::VectortoPoint(const VectorC &vec)
{
	Array1dC<Index2dC> out((SizeT)(vec.Size()/2));
	for(SizeT i =0; i < (vec.Size()/2); i++)
	{
		Index2dC inp(vec[i +(SizeT)(vec.Size()/2)],vec[i]);
		out[i] = inp.Copy();
	}
	return out;
}

Array1dC<Point2dC> JMeasurementModel::VectortoPointArray(const VectorC &vec)
{
	Array1dC<Point2dC> out((SizeT)(vec.Size()/2));
	for(SizeT i =0; i < (vec.Size()/2); i++)
	{
		Point2dC inp(vec[i +(SizeT)(vec.Size()/2)],vec[i]);
		out[i] = inp.Copy();
	}
	return out;
}

ImageC<RealT> JMeasurementModel::ComputeJImage(const ImageC<RealRGBValueC> &img, const MeanCovarianceC &mc)
{
	RealT logsigma = Log((RealT)mc.Covariance().Det());
	ImageC<RealT> out(img.Frame());
	out.Fill(0.0);
	for(Array2dIter2C<RealRGBValueC, RealT> it(img, out); it; it++)
	{
		RealRGBValueC rgb = it.Data1();
		VectorC pt(rgb.Red(), rgb.Green(), rgb.Blue());
		it.Data2() = -logsigma - mc.MahalanobisDistance(pt);
	}
	return out;
}
ImageC<RealT> JMeasurementModel::ComputeJImage(const ImageC<VectorC> &img, const MeanCovarianceC &mc)
{
	RealT logsigma = Log((RealT)mc.Covariance().Det());
	ImageC<RealT> out(img.Frame());
	out.Fill(0.0);
	for(Array2dIter2C<VectorC, RealT> it(img, out); it; it++)
	{
		VectorC rgb = it.Data1();
		it.Data2() = -logsigma - mc.MahalanobisDistance(rgb);
	}
	return out;
}
/*
RealT JMeasurementModel::GetJWeight(const Array1dC<Index2dC> &ar, ImageC<RealT> &mc, ImageC<RealT> &sk)
{
	Index2dC ctr = ar.Sum().Copy();
	ctr.Row() /= ar.Size();
	ctr.Col() /= ar.Size();
	RealT sumsquare = 0.0;
	for(Array1dIterC<Index2dC> it(ar); it; it++)
	{
		LinePP2dC line(ctr,(*it));
		Array1dC<RealT> jarray(window);
		jarray.Fill(0.0);
		IndexC arind = 0;
		for(RealT i = 0.0; i < (window/10.0); )
		{
			Point2dC loc = line.Point(i);
			if(jlip.Frame().Contains(loc))
			{
				jarray[arind] = jlip[loc] - jskin[loc]; 
			}
			i += 0.1;
			arind ++;
		}
		cout<<"JARRAY "<<jarray<<endl;
		//now we have an array of window measurements, starting from the outside, look for the first zero crossing
		arind = 0;
		IndexC ind = 10;
		for(IndexC i = arind; i < window; i++)
		{
			if((jarray[i] > 0.0 )&&(jarray[i + 1] < 0.0))
			{
				ind = i;
				cout<<"Index Chosen = "<<i<<endl;
			}
		}
		Index2dC lipbndry = (Index2dC)line.Point((RealT)(ind/10.0));
		sumsquare += (RealT)(*it).SqrEuclidDistance(lipbndry);
	}
	if(sumsquare == 0.0)
		return 1.0;
	else 
		return (1.0/sumsquare);
}

RealT JMeasurementModel::GetJWeight(const Array1dC<Index2dC> &ar, ImageC<RealT> &mc, ImageC<RealT> &sk)
{
	Index2dC ctr = ar.Sum().Copy();
	ctr.Row() /= ar.Size();
	ctr.Col() /= ar.Size();
	RealT sumsquare = 0.0;
	for(Array1dIterC<Index2dC> it(ar); it; it++)
	{
		//Can use vertical jdiff image
		LinePP2dC line(ctr,(*it));
		Array1dC<RealT> jarray(window);
		jarray.Fill(0.0);
		IndexC arind = 0;
		for(RealT i = Abs(1.0 -(window/20.0)) ; i < (1+(window/20.0)); )//(Arr start - half od window BEFORE 1.0, hence div by 20)
		{
			Point2dC loc = line.Point(i);
			if(jlip.Frame().Contains(loc))
			{
				jarray[arind] = vertical[loc]; //Insert vertical gradient 
			}
			i += 0.1;
			arind ++;
		}
		//cout<<"JARRAY "<<jarray<<endl;
		//now we have an array of window measurements, look for largest positive value
		arind = 0;
		IndexC ind = 10;
		for(IndexC i = arind; i < window-1; i++)
		{
			RealT absval = Abs(jarray[i]);
			cout<<"ABS VAL = "<<absval<<endl;
			if((jarray[i] > 0)&&(jarray[i+1] < 0))
			{
				ind = i;
				//val = absval;
				//cout<<"Index Chosen = "<<i<<"New Value = "<<val<<endl;
			}
		}
		Index2dC lipbndry = (Index2dC)line.Point((RealT)(ind/10.0));
		sumsquare += (RealT)(*it).SqrEuclidDistance(lipbndry);
	}
	if(sumsquare == 0.0)
		return 1.0;
	else 
		return (1.0/sumsquare);
}*/

RealT JMeasurementModel::GetJWeight(const Array1dC<Point2dC> &ar, ImageC<RealT> &mc, ImageC<RealT> &sk)
{
	//Using Points, Get B-Spline Control Points
	UIntT order = 3; UIntT ncp = 11;
	BSplineC bspl(order,ncp, BSplineC::UOPEN);
	bspl.SetControlPoints(ar);
	//Now Measure the distance
	ImageC<RealT> j_im = mc - sk;
	UIntT curve_res = 20;
	RealT incr = 1.0/(RealT)curve_res;	
	RealT j_wt = 0.0;IntT search = window;
	for(RealT pval = 0; pval < curve_res; pval ++)		
	{
		RealT par = (RealT)pval*incr;
		LinePP2dC normal = bspl.GetNormalLine(par); //Get B-Spline Normal Line
		Point2dC start = normal.Point(0.0);
		//Now just search from the start to the end for the binary condition to end
		RealT index;
		for(RealT i = -5; i < search; i++)
		{
			Point2dC cur = normal.Point(i); Point2dC next = normal.Point(i+1);
			if(j_im.Frame().Contains(next))
			{
				if((j_im[cur] >= 0)&&(j_im[next] < 0))
				{
					index = i;
				}
			}
		}
		j_wt += GetGaussianValue(3,index);
	}
	return j_wt;
}


RealT JMeasurementModel::GetJWeight(const Array1dC<Point2dC> &ar)
{
	ImageC<VectorC> orig_vect = GetNormRGVals(origimg);
	//Using Points, Get B-Spline Control Points
	UIntT order = 3; UIntT ncp = 11;
	BSplineC bspl(order,ncp, BSplineC::UOPEN);
	bspl.SetControlPoints(ar);
	//Now Measure the distance
	UIntT curve_res = 20;
	RealT incr = 1.0/(RealT)curve_res;	
	RealT j_wt = 0.0;IntT search = window;
	for(RealT pval = 0; pval < curve_res; pval ++)		
	{
		RealT par = (RealT)pval*incr;
		LinePP2dC normal = bspl.GetNormalLine(par); //Get B-Spline Normal Line
		Point2dC start = normal.Point(0.0);
		//Now just search from the start to the end for the binary condition to end
		//Change condition to start looking outwards in
		RealT index;
		//for(RealT i = -5; i < search; i++)
		for(RealT i = search; i > 0; i--)
		{
			//Point2dC cur = normal.Point(i); Point2dC next = normal.Point(i+1);
			Point2dC cur = normal.Point(i); Point2dC next = normal.Point(i-1);
			if(orig_vect.Frame().Contains(next))
			{
				RealT j_cur = ComputeJ(lip,orig_vect[cur]) - ComputeJ(skin,orig_vect[cur]);
				RealT j_next = ComputeJ(lip,orig_vect[next]) - ComputeJ(skin,orig_vect[next]);
				//if((j_cur >= 0)&&(j_next < 0))
				if((j_cur <= 0)&&(j_next > 0))
				{
					index = i;
					break;
				}
			}
		}
		j_wt += GetGaussianValue(3,index);
	}
	return j_wt;
}
RealT JMeasurementModel::ComputeJ(const MeanCovarianceC &mc, const VectorC &vec)
{
	return (-Log(mc.Covariance().Det()) - mc.MahalanobisDistance(vec));
}
RealT JMeasurementModel::GetGaussianValue(const RealT &sigma, const RealT &xval)
{
	//Generate a 1D gaussian y-value given the value of (x - mu)
	return (1/(sigma*Sqrt(2*pi)))*Exp(-(Sqr(xval)/(2*Sqr(sigma))));
}

