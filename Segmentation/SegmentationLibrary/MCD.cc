//File: MCD.cc
//Class - defines the interface to the Minimum Covariance Determinant Function
// Takes in as an input, a set of Sample of Vectors. Performs MCD estimation on that set of vectors.
// Can return various final estimates from procedure. 
//Author: Bud Goswami
//Date: 29.01.09

//Datavariables:
//Private Data Members
//Algorithm Parameters
//RealT hvalue; //to store confidence in parameters
//UIntT nvalue; //to store number of iterations to convergence
//Starting data population: DListC<Tuple2C<VectorC,Index2dC> > start_pop;
//Elemental Set: Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > elemental_set;
//Final Set: Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > end_pop;


#include "MCD.hh"

//ImageC<RealRGBValueC> related methods - start with conversion and then define constructors
DListC<Tuple2C<VectorC,Index2dC> > MCD::GetStartPopFromImage(const ImageC<RealRGBValueC> &img)
{
	DListC<Tuple2C<VectorC,Index2dC> > res;
	Tuple2C<VectorC, Index2dC> data;
	for(Array2dIterC<RealRGBValueC> it(img); it; it++)
	{
		data.Data1() = VectorC((*it).Red(),(*it).Green(),(*it).Blue()).Copy();
		data.Data2() = it.Index();
		res.Append(data.Copy());
	}
	return res;
}

MCD::MCD(const ImageC<RealRGBValueC> &rgbimg, const RealT &hval, const UIntT &numiter)
{
	start_pop = GetStartPopFromImage(rgbimg);
	if(hval >= 0.5)
		hvalue = hval;
	else
		cerr<<"Breakdown point for this estimator is 0.5, incorrect h value entered"<<endl;
	nvalue = numiter; 	
}
MCD::MCD(const ImageC<RealRGBValueC> &rgbimg, const RealT &hval)
{
	start_pop = GetStartPopFromImage(rgbimg);
	if(hval >= 0.5)
		hvalue = hval;
	else
		cerr<<"Breakdown point for this estimator is 0.5, incorrect h value entered"<<endl;
	nvalue = 10;
}
MCD::MCD(const ImageC<RealRGBValueC> &rgbimg, const UIntT &numiter)
{
	start_pop = GetStartPopFromImage(rgbimg);
	nvalue = numiter;
	hvalue = 0.75;
}
MCD::MCD(const ImageC<RealRGBValueC> &rgbimg)
{
	start_pop = GetStartPopFromImage(rgbimg);
	nvalue = 10;
	hvalue = 0.75;
}

//ImageC<VectorC> related methods - in case you want to instantiate with an image of data-type Vectors
DListC<Tuple2C<VectorC,Index2dC> > MCD::GetStartPopFromImageVector(const ImageC<VectorC> &img)
{
	DListC<Tuple2C<VectorC,Index2dC> > res;
	Tuple2C<VectorC, Index2dC> data;
	for(Array2dIterC<VectorC> it(img); it; it++)
	{
		data.Data1() = (*it).Copy();
		data.Data2() = it.Index();
		res.Append(data.Copy());
	}
	return res;
}
MCD::MCD(const ImageC<VectorC> &vecimg, const RealT &hval, const UIntT &numiter)
{
	start_pop = GetStartPopFromImageVector(vecimg);
	if(hval >= 0.5)
		hvalue = hval;
	else
		cerr<<"Breakdown point for this estimator is 0.5, incorrect h value entered"<<endl;
	nvalue = numiter; 	
}
MCD::MCD(const ImageC<VectorC> &vecimg, const RealT &hval)
{
	start_pop = GetStartPopFromImageVector(vecimg);
	if(hval >= 0.5)
		hvalue = hval;
	else
		cerr<<"Breakdown point for this estimator is 0.5, incorrect h value entered"<<endl;
	nvalue = 10; 	

}
MCD::MCD(const ImageC<VectorC> &vecimg, const UIntT &numiter)
{
	start_pop = GetStartPopFromImageVector(vecimg);
	hvalue = 0.75;
	nvalue = numiter;
}
MCD::MCD(const ImageC<VectorC> &vecimg)
{
	start_pop = GetStartPopFromImageVector(vecimg);
	hvalue = 0.75;
	nvalue = 10;	
}
//Constructors with SampleC<VectorC>, RealT and UIntT
//NOTE that using SampleDataC for images assumes that the vector includes index information
//Form: vec[0] = imgdata0...vec[dim-1] = imgdata_dim; vec[dim] = Index.Row(), vec[dim+1] = Index.Col() 
DListC<Tuple2C<VectorC,Index2dC> > MCD::GetStartPopFromSampleData(SampleData &img)
{
	DListC<Tuple2C<VectorC,Index2dC> > res; 
	DListC<VectorC> orig = img.GetSmple().Copy();
	IntT dims = img.GetDimensions();
	for(DLIterC<VectorC> it(orig); it; it++)
	{
		VectorC data(dims);
		for(IntT i = 0; i < dims; i++)
		{
			data[i] = (*it)[i];
		}
		Index2dC ind((*it)[dims],(*it)[dims+1]);
		Tuple2C<VectorC, Index2dC> res_elem(data.Copy(), ind.Copy());
		res.Append(res_elem.Copy());
	}
	return res;
}
MCD::MCD(SampleData &smpl, const RealT &hval, const UIntT &numiter)
{
	start_pop = GetStartPopFromSampleData(smpl);
	if(hval >= 0.5)
		hvalue = hval;
	else
		cerr<<"Breakdown point for this estimator is 0.5, incorrect h value entered"<<endl;
	nvalue = numiter; 	
}
MCD::MCD(SampleData &smpl, const RealT &hval)
{
	start_pop = GetStartPopFromSampleData(smpl);
	if(hval >= 0.5)
		hvalue = hval;
	else
		cerr<<"Breakdown point for this estimator is 0.5, incorrect h value entered"<<endl;
	nvalue = 10;	
}
MCD::MCD(SampleData &smpl, const UIntT &numiter)
{
	start_pop = GetStartPopFromSampleData(smpl);
	hvalue = 0.75;
	nvalue = numiter; 		
}
MCD::MCD(SampleData &smpl)
{
	start_pop = GetStartPopFromSampleData(smpl);
	hvalue = 0.75;
	nvalue = 10; 		
}

//Constructors for any 3-dimensional pixel type in RAVL (As long as they are "Real")
//E.g. RealHSVValueC, RealYUVValueC...
DListC<Tuple2C<VectorC,Index2dC> > MCD::GetStartPopFromImageVector(const ImageC<TFVectorC<RealT,3> > &img)
{
	DListC<Tuple2C<VectorC,Index2dC> > res;
	Tuple2C<VectorC, Index2dC> data;
	for(Array2dIterC<TFVectorC<RealT,3> > it(img); it; it++)
	{
		VectorC vec(3);
		for(UIntT i = 0; i < 3; i++)
		{
			vec[i] = (*it)[i];
		}
		data.Data1() = vec.Copy();
		data.Data2() = it.Index();
		res.Append(data.Copy());
	}
	return res;	
}

MCD::MCD(const ImageC<TFVectorC<RealT,3> > &rgbimg, const RealT &hval, const UIntT &numiter)
{
	start_pop = GetStartPopFromImageVector(rgbimg);
	if(hval >= 0.5)
		hvalue = hval;
	else
		cerr<<"Breakdown point for this estimator is 0.5, incorrect h value entered"<<endl;
	nvalue = numiter;
}
	
MCD::MCD(const ImageC<TFVectorC<RealT,3> > &rgbimg, const RealT &hval)
{
	start_pop = GetStartPopFromImageVector(rgbimg);
	if(hval >= 0.5)
		hvalue = hval;
	else
		cerr<<"Breakdown point for this estimator is 0.5, incorrect h value entered"<<endl;
	nvalue = 10; 		
}

MCD::MCD(const ImageC<TFVectorC<RealT,3> > &rgbimg, const UIntT &numiter)
{
	start_pop = GetStartPopFromImageVector(rgbimg);
	hvalue = 0.75;
	nvalue = numiter;
}
	
MCD::MCD(const ImageC<TFVectorC<RealT,3> > &rgbimg)
{
	start_pop = GetStartPopFromImageVector(rgbimg);
	hvalue = 0.75;
	nvalue = 10; 	
}
	
//GetElementalSet()
//This generates a uniformly random number and assigns it to the Elemental Set
//If the covariance of the elemental set (which must be of dimensionality one greater than data vector) is
//zero, then elements are added to this set until the covariance is non-zero
Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > MCD::GetElementalSet(void)
{
	UIntT p = start_pop.First().Data1().Size(); //to store the dimensionality of the data, p
	UIntT totelem = start_pop.Size(); //to store the population of data, n
	DListC<Tuple2C<VectorC,Index2dC> > elem_set; // to store the actual elemental set
   SumsNd2C elem_data(p); //to store the data in the elemental set
	bool inclusive = false; //for the random number inclusion
	//Draw p+1 uniformly random observations from data to form elemental set
	for(UIntT i = 0; i < p+1; i++)
	{
		UIntT ind = (UIntT)(Random1(inclusive)*totelem);
		elem_set.Append(start_pop.Nth(ind).Copy());
		elem_data += elem_set.Last().Data1();
	}
	//Compute covariance and determinant
	bool sampleStatistics = false;
	MeanCovarianceC elem_mc = elem_data.MeanCovariance(sampleStatistics);
	//If covariance determinant is zero, add elements until it is non-zero
	while(elem_mc.Covariance().Det() == 0.0)
	{
		UIntT ind = (UIntT)(Random1(inclusive)*totelem);
		elem_set.Append(start_pop.Nth(ind).Copy());
		elem_data += elem_set.Last().Data1();
		elem_mc = elem_data.MeanCovariance(); 
	}
	Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > out(elem_mc, elem_set);
	return out;
}
//mscomp function for MergeSort function for the DList generated in the GetCStep function below
bool mscomp(const Tuple3C<VectorC,Index2dC, RealT> &dat1, const Tuple3C<VectorC,Index2dC, RealT> &dat2)
{
	return (dat1.Data3() <= dat2.Data3() ? true:false);
}
//GetCStep()
//Performs MCD CStep estimation as described by URS in Ch5
//Given initial MC estimate
//Compute Mahalanobis distances of all data points from the data
//Merge-sort the data by Mahalanobis distances and get new data set
//Use new data set to obtain updated value of MeanCovariance
MeanCovarianceC MCD::GetCStep(const MeanCovarianceC &mcstart, const RealT &h, const  DListC<Tuple2C<VectorC,Index2dC> > &st_pop)	
{
	DListC<Tuple3C<VectorC,Index2dC, RealT> > weighted_pop;
	DListC<Tuple2C<VectorC,Index2dC> > end_pop_data;
	//For all the elements in the start_population, compute the Mahalanobis distance from estimate
	for(DLIterC<Tuple2C<VectorC,Index2dC> > it(st_pop); it; it++)
	{
		RealT ms_dist = mcstart.MahalanobisDistance((*it).Data1());
		Tuple3C<VectorC,Index2dC, RealT> elem((*it).Data1(),(*it).Data2(),ms_dist);
		weighted_pop.Append(elem.Copy());
	}
	//cout<<"Weighted Pop - "<<weighted_pop.Size()<<"\t First Elem - "<<weighted_pop.First()<<endl;
	//Now, given a data set with computed Mahalanobis Distances, merge sort it for size
	weighted_pop.MergeSort(mscomp);
	//Parse through the first h fraction of the data set and use it for refined meancov estimation
	SumsNd2C elem_data(st_pop.First().Data1().Size()); //This will get the result of the MeanCov computation
	UIntT i = 0;
	for(DLIterC<Tuple3C<VectorC,Index2dC,RealT> > it(weighted_pop); it; it++)
	{
		if(i > (UIntT)(h*weighted_pop.Size()))
			break;
		else
		{
			elem_data += (*it).Data1();
			Tuple2C<VectorC,Index2dC> end_pop_elem((*it).Data1(),(*it).Data2());
			end_pop_data.Append(end_pop_elem.Copy()); //insert element into end_pop_data structure
			i++;
		}
	}
	bool SampleStats = false; //Used to compute MeanCovariance
	MeanCovarianceC mcref_est = elem_data.MeanCovariance(SampleStats); //Compute MeanCovarianceC
	//cout<<"Computed Refinement = "<<mcref_est<<endl;
	//Assign results to final_pop and also return the estimate
	end_pop.Data1() = mcref_est.Copy();
	end_pop.Data2() = end_pop_data.Copy();
	return mcref_est;
}

//Apply Function
//Performs the MCD estimation as stated out in URS thesis Ch 5
//Using start_pop, obtain Elemental Set and Start Estimates of MeanCovariance
//Using this start estimate, iteratively apply the GetCStep() function until convergence
//Convergence conditions - MeanCovariance.Covariance.Det stops changing or equal to zero
//Alternative convegence condition (for tractability) - after numiter iterations (usually=10)
//Eccentricity in Method*********
//For convergence condition, we only want the equality to be to three s.f.
//In number representations, this is not possible. Therefore, for comparison condition:
//v1*1000 == v2*1000 
Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > MCD::Apply(void)
{
	//P.S. In Single Pass MCD, important NOT to change private data members
	//For multi-pass MCD (e.g. Cascaded/ K-S MCD, derive classes and also implement Set() accessor functions
	//Additionally, can implement new Apply Function (assuming they are children and not aggregate classes with MCD as component)
	elemental_set = GetElementalSet();
	//cout<<"Elemental Set Obtained\t"<<elemental_set.Data1()<<endl;
   MeanCovarianceC start_ref_est(elemental_set.Data1().Copy());
	//Run C-Step nvalue number of times
	for(UIntT i = 0; i <= nvalue; i++)
	{
		MeanCovarianceC fin_ref_est;
		fin_ref_est = GetCStep(start_ref_est, hvalue, start_pop);
		//cout<<"Inside C_Step loop - "<<i<<"\n fin_ref_est.Covariance().Det() = "<<fin_ref_est.Covariance().Det()<<endl;;
		//<<"\n start_ref_est.Covariance().Det()="<<start_ref_est.Covariance().Det()<<endl;
		
		//Check for convergence (READ ECCENTRICITY NOTE IN METHOD)
		//REASON FOR MULTIPLYING BY 1000 IS TO COMPARE TO 3 S.F.
		if(((UIntT)(fin_ref_est.Covariance().Det()*100.0) == (UIntT)(start_ref_est.Covariance().Det()*100.0))||(fin_ref_est.Covariance().Det() == 0.0))
			break;
		else
			start_ref_est = fin_ref_est.Copy();
	}
	return end_pop;
}
