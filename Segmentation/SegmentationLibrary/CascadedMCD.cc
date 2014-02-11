//File: CascadedMCD.cc
//Class - defines the interface to the CascadedMinimum Covariance Determinant Function
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


#include "CascadedMCD.hh"

//CascadedMCD Apply Function
//Performs the MCD estimation as stated out in URS thesis Ch 5 iteratively
//Using start_pop, obtain Elemental Set and Start Estimates of MeanCovariance
//Using this start estimate, iteratively apply the GetCStep() function until convergence at each step
//After convergence, remove the final pop from the data and reapply
//Convergence conditions - MeanCovariance.Covariance.Det stops changing or equal to zero
//Eccentricity in Method*********
DListC< Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > CascadedMCD::Apply(void)
{
	DListC< Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > output;
	DListC<Tuple2C<VectorC,Index2dC> > current_start = start_pop.Copy();
	do
	{
		//cout<<"Start_pop size = "<<current_start.Size()<<endl;;
		Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > cur_end = MCD::Apply();
		output.Append(cur_end.Copy());
		//cout<<"End_pop Size = "<<cur_end.Data2().Size()<<endl;
		//now end_pop will have initial end_population
		//remove this population from current_start
		for(DLIterC<Tuple2C<VectorC,Index2dC> > it(cur_end.Data2()); it; it++)
		{
			current_start.Del(*it);
		}
		//cout<<"New start_pop size = "<<current_start.Size()<<endl; 
		
		//set start_pop to be equal to current_start
		SetStartPop(current_start.Copy());
	}while(current_start.Size() >= minsamp);
	return output;
}
