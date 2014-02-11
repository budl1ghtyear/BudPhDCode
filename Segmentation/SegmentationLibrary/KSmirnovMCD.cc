//File: KSmirnovMCD.cc
//Class - defines the interface to the KSmirnovMinimum Covariance Determinant Function
// Takes in as an input, a set of Sample of Vectors. Performs iterative MCD estimation on that set of vectors until minimum sample size is reached.
// Can return various final estimates from procedure. 
//Author: Bud Goswami
//Date: 16.02.09

//Datavariables:
//Private Data Members
//Algorithm Parameters
//RealT hvalue; //to store confidence in parameters
//UIntT nvalue; //to store number of iterations to convergence
//UIntT minsamp; // to store the minimum number of samples with which to quit the algorithm
//Starting data population: DListC<Tuple2C<VectorC,Index2dC> > start_pop;
//Elemental Set: Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > elemental_set;
//Final Set: Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > end_pop;


#include "KSmirnovMCD.hh"
//GetZerothEstimate()
//Given an estimate using MCD, calculate refined "zeroth" estimates from these
//T remains unchanged, S = reweighted S
//Reqeight - maxmsdist/chi squared correction factor
//Given a fin_pop estimate, this will obtain the T_0 and S_0 estimates from the data
MeanCovarianceC KSmirnovMCD::GetZerothEstimate(const MeanCovarianceC &mc, const RealT &h, const RealT &dh, const UIntT &mval)
{
	MeanCovarianceC out = mc.Copy();
	out.Covariance() *= (dh/(RealT)chisquarecdistribution(mval,h));
	return out;
}

//FinMCEstimate()
//Given the original MCD MC estimate and the Zeroth refined MCD estimate, reweight the observations
//Now compute Mean and Covariance using SumsN2dC
MeanCovarianceC KSmirnovMCD::FinMCEstimate(const MeanCovarianceC &mc, const MeanCovarianceC &mczero, const UIntT &m, const RealT &h)
{
	SumsNd2C mccal(m);
	bool sampleStatistics = true;
	RealT chival = (RealT)chisquarecdistribution(m,h);
	//DListC<Tuple2C<VectorC,Index2dC> > start_pop; - iterate through this and calculate weights
	for(DLIterC<Tuple2C<VectorC,Index2dC> > it(start_pop); it; it++)
	{
		RealT w = mczero.MahalanobisDistance((*it).Data1());	
		if(w <= chival)
			w = 1.0;
		else
			w = 0.0;
		RealT q = mc.MahalanobisDistance((*it).Data1());
		mccal.Add(((*it).Data1()),(RealT)w*q);
	}
	return mccal.MeanCovariance(sampleStatistics);
}
//Given a specific set of measurements, compute the K-S Test
//First, compute an array of intergral sums for Fobs and also estimate i_max given end_pop
RealT KSmirnovMCD::KSTest(const MeanCovarianceC &mc, const UIntT &m)
{
	//Create two arrays of F values
	Array1dC<RealT> fobs(end_pop.Data2().Size()); fobs.Fill(0.0);
	Array1dC<RealT> fthr(end_pop.Data2().Size()); fobs.Fill(0.0);
	//End_pop is sorted in increasing order, so we can just use the observations in the order in which they arrive
	//Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > end_pop;
	//First pass, fill out Fthr and Num of Fobs
	IntT index = 0;
	RealT num = 0.0;
	for(DLIterC<Tuple2C<VectorC,Index2dC> > it(end_pop.Data2()); it; it++)
	{
		RealT qj = mc.MahalanobisDistance((*it).Data1());
		RealT chival = (RealT)chisquarecdistribution(m,0.975);
		//cout<<"Computed Chi Square Distrib - "<<chival<<"\t qj = "<<qj<<endl;
		if(index == 0)
		{	
			fobs[0] = 0.0;
			fthr[0] = 0.0;
		}
		else
		{
			if(qj < chival)
			{
				num += qj; //create a running sum and insert it as fobs value
				fobs[index] = num;
				fthr[index] = (RealT)invchisquaredistribution(m,qj);
			}	
		}
		index++;
	}
	fobs /= num;
	RealT ksval = 0.0;
	for(Array1dIter2C<RealT, RealT> it(fobs, fthr); it; it++)
	{
		RealT dif = Abs(it.Data1() - it.Data2());
		if(dif > ksval)
			ksval = dif;
	}
	return ksval;
}
//KSmirnovMCD Apply Function
//Performs the MCD estimation as stated out in URS thesis Ch 5 iteratively
//Using start_pop, obtain Elemental Set and Start Estimates of MeanCovariance
//Using this start estimate, iteratively apply the GetCStep() function until convergence at each step
//After convergence, remove the final pop from the data and reapply
//Convergence conditions - MeanCovariance.Covariance.Det stops changing or equal to zero
//Eccentricity in Method*********
DListC< Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > KSmirnovMCD::Apply(void)
{
	DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > output;
	//Create a for loop to vary the hvalue between 0.6 and 0.9 in increments of 0.1
	//cout<<"Inside Apply Method"<<endl;
	UIntT dims = start_pop.First().Data1().Size();
	DListC<Tuple2C<VectorC,Index2dC> > cur_start_pop = start_pop.Copy();
	Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > dh_end_pop; //to store the best end_pop
	do
	{
		//cout<<"Inside Apply Method - starting new cluster"<<endl;
		//New cluster estimation iteration starts here
		RealT Dhvalue = 0.0; //to check for the best cluster
		RealT h = 0.6; //alternate hypothesis for h
		while(h < 1.0)
		{
			//cout<<"\nH Value - "<<h<<"\t Data Population - "<<start_pop.Size()<<"\t";
			//Set the MCD Value of h
			SetH(h);
			MCD::Apply();
			//cout<<"Data Size of End_pop = "<<end_pop.Data2().Size()<<endl;
			MeanCovarianceC mcest = end_pop.Data1().Copy();
			//cout<<"Mean Covariance Init - "<<mcest.Covariance().Det()<<endl;
			RealT dh = mcest.MahalanobisDistance(end_pop.Data2().Last().Data1());
			MeanCovarianceC mczero = GetZerothEstimate(mcest,h, dh, dims); //get zeroth refinment
			//cout<<"Mean Covariance Zero - "<<mczero.Covariance().Det()<<endl;
			MeanCovarianceC mcfin = FinMCEstimate(mcest,mczero,dims,h); //compute least squares MC vals
			//cout<<"Mean Covariance Final - "<<mcfin.Covariance().Det()<<endl;
			//Now, given new estimate, insert the refined estimate and obtained KS test result
			RealT ksval = KSTest(mcfin,dims);
			//cout<<"Current Value for Dh = "<<Dhvalue<<"\tKS Test Result = "<<ksval<<"\t H = "<<h<<endl;
			if (ksval >= Dhvalue)
			{
				Dhvalue = ksval; 
				dh_end_pop = end_pop.Copy();
			}
			h += 0.1;
		}
		//Copy best end_population to the output result
		output.Append(dh_end_pop.Copy());
		//At this point, we have a set of data - dh_end_pop which we can use to eliminate the start_pop observations
		for(DLIterC<Tuple2C<VectorC,Index2dC> > it(dh_end_pop.Data2()); it; it++)
		{
			cur_start_pop.Del(*it);
		}
		start_pop = cur_start_pop.Copy();
	}while(cur_start_pop.Size() > minsamp);
	//cout<<"Exhausted Data"<<endl;
	return output;
}
