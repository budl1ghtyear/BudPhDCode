#include"SpectralClustering.hh"

//SpectralClustering.cc
//Declare the functionality of the spectral clustering class
//Author - Bud Goswami
//Date - 17.02.09

//Functions for constructor usage
Tuple2C<SampleC<VectorC>, DListC<Tuple2C<VectorC, Index2dC> > > SpectralClustering::ConvertImageToSamples(const ImageC<RealRGBValueC> &rgb)
{
	SArray1dC<VectorC> d(rgb.Frame().Area());
	DListC<Tuple2C<VectorC, Index2dC> > list;
	IntT ind = 0;
	for(Array2dIterC<RealRGBValueC> it(rgb); it; it++)
	{
		VectorC vec((*it).Red(),(*it).Blue(),(*it).Blue());
		d[ind++] = vec.Copy();
		list.Append(Tuple2C<VectorC, Index2dC>(vec.Copy(), it.Index()).Copy());
 	}
	SampleC<VectorC> sample(d.Copy());	
	Tuple2C<SampleC<VectorC>, DListC<Tuple2C<VectorC, Index2dC> > > out(sample, list);
	return out;
}


Tuple2C<SampleC<VectorC>, DListC<Tuple2C<VectorC, Index2dC> > > SpectralClustering::ConvertImageToSamples(const ImageC<VectorC> &rgb)
{
	SArray1dC<VectorC> d(rgb.Frame().Area());
	DListC<Tuple2C<VectorC, Index2dC> > list;
	IntT ind = 0;
	for(Array2dIterC<VectorC> it(rgb); it; it++)
	{
		VectorC vec((*it).Copy());
		d[ind++] = vec.Copy();
		list.Append(Tuple2C<VectorC, Index2dC>(vec.Copy(), it.Index()).Copy());
 	}
	SampleC<VectorC> sample(d.Copy());	
	Tuple2C<SampleC<VectorC>, DListC<Tuple2C<VectorC, Index2dC> > > out(sample, list);
	return out;
}

SpectralClustering::SpectralClustering(const ImageC<RealRGBValueC> &rgbimg, const UIntT &cl)
{
	imrect = rgbimg.Frame();
	Tuple2C<SampleC<VectorC>, DListC<Tuple2C<VectorC, Index2dC> > > out = ConvertImageToSamples(rgbimg);
	obs = out.Data1().Copy();
	obsind = out.Data2().Copy();
	numclust = cl;
}
SpectralClustering::SpectralClustering(const ImageC<RealRGBValueC> &rgbimg)
{
	imrect = rgbimg.Frame();
	Tuple2C<SampleC<VectorC>, DListC<Tuple2C<VectorC, Index2dC> > > out = ConvertImageToSamples(rgbimg);
	obs = out.Data1().Copy();
	obsind = out.Data2().Copy();
	numclust = 4; //Default number of clusters (remember to add one for 1 eigenvector)
}
SpectralClustering::SpectralClustering(const ImageC<VectorC> &vecimg, const UIntT &cl)
{
	imrect = vecimg.Frame();
	Tuple2C<SampleC<VectorC>, DListC<Tuple2C<VectorC, Index2dC> > > out = ConvertImageToSamples(vecimg);
	obs = out.Data1().Copy();
	obsind = out.Data2().Copy();
	numclust = cl; //Default number of clusters (remember to add one for 1 eigenvector)
}
SpectralClustering::SpectralClustering(const ImageC<VectorC> &vecimg)
{
	imrect = vecimg.Frame();
	Tuple2C<SampleC<VectorC>, DListC<Tuple2C<VectorC, Index2dC> > > out = ConvertImageToSamples(vecimg);
	obs = out.Data1().Copy();
	obsind = out.Data2().Copy();
	numclust = 4; //Default number of clusters (remember to add one for 1 eigenvector)
}

//Compute Similarity Connected Graph Depending on Similarity Function
Tuple2C<MatrixC,MatrixC> SpectralClustering::ComputeSDMatrix(SimType s)
{
	MatrixC sim(obs.Size(), obs.Size(), 0.0);
	MatrixC diag(obs.Size(), obs.Size(), 0.0);
	switch(s)
	{
		case GAUSSIANKERNEL:
		{
			//cout<<"Computing Gaussian Similarity & Diagonal Transform"<<endl;
			for(UIntT i = 0; i < obs.Size(); i++)
			{	
				RealT rowsum = 0.0;
				VectorC xi = obs[i].Copy();
				for(UIntT j = 0; j < obs.Size(); j++)
				{
					if(j != i)
					{
						//if(j >= i)
						{
							VectorC xj = obs[j].Copy();
							//Compute Gaussian Kernel Value with sigma 0.5
							//Formula is -||x_i - x_j||^2 / 2sigma^2 
							//Perform computation for Similarity Matrix
							sim[i][j] = Exp(-(xi.SqrEuclidDistance(xj)/0.5));
							//sim[j][i] = sim[i][j];
						}
					}
					rowsum += sim[i][j];
				}
				//Perform computation for Diagonal Row-Sum Matrix
				diag[i][i] = rowsum;
			}
		}
		break;
		
		case SHINCUT:
		{
			for(UIntT i = 0; i < obs.Size(); i++)
			{
				VectorC fi = obs[i];
				Index2dC xi = obsind.Nth(i).Data2();
				RealT rowsum = 0.0;
				for(UIntT j = 0; j < obs.Size(); j++)
				{	
					Index2dC xj = obsind.Nth(j).Data2();
					RealT eucx = xi.Distance(xj,SQUARE_EUCLID);
					//r = 10
					if(eucx > 10)
					{
						//sigma_x = 4.0, sigma_I = 0.01
						eucx = Exp(eucx * -0.25);					
						RealT fin = eucx * Exp(fi.SqrEuclidDistance(obs[j]) * -100.0);
						sim[i][j] = fin;						
					}
					rowsum += sim[i][j];
				}
				diag[i][i] = rowsum;		
			}	
		}
		break;
	}
	//cout<<"Finished Computing Gaussian Similarity & Diagonal Transform"<<endl;
	return Tuple2C<MatrixC, MatrixC> (sim,diag);
}

MatrixC SpectralClustering::ComputeLaplacian(const Tuple2C<MatrixC, MatrixC> &sd, LapType L)
{
	//cout<<"Computing Laplacian Matrix"<<endl;
	MatrixC lmat((MatrixC)sd.Data1());
	//Smart stuff to compute quick subtraction by a diagonal matrix
	//Only diagonal elements are subtracted
	//Remaining part of matrix and just 
	lmat *= (-1.0);
	for(UIntT i = 0; i < sd.Data1().Rows(); i ++)
	{
		lmat[i][i] +=  sd.Data2()[i][i];
	}
	//At this point, L = D - S;
	//Now compute the laplacian depending on what type of normalisation to use
	switch(L)
	{
		case RANDOMWALK:
		{
			lmat = sd.Data1().Inverse()*lmat;
		}
		break;
		
		case UNNORMALISED:
		{
			
		}
		break;
		
		case SYMMETRIC:
		{
				//Here we need to implement - L = D^(-1/2)SD^(-1/2)
				cerr<<"Method Not Implemented Yet, Use Standard Laplacian"<<endl;
		}
		break;
	}
	//cout<<"Finished Computing Laplacian Matrix"<<endl;
	cout<<lmat<<endl;
	return lmat;	
}

MatrixC SpectralClustering::GetNormalisedEigenVectors(MatrixC &lmat)
{
	/* //This method uses standard eigvector call from RAVL
	MatrixC eig; //to store the eigenvectors
	VectorC eigval = EigenVectors(lmat,eig); //to store the eigenvalues
	cout<<"Computed Eigenvalues and EigenVectors "<<endl;
	eig.NormaliseRows(); //now we have a normalised row set of vertical eigenvectors
	cout<<"Computed Normalised Rows for EigenVectors"<<endl;
	MatrixC out = eig.SubMatrix(lmat.Rows(), numclust);
	*/
	
	//Using fast eigenvector call
	//VectorC testval = FastEigenValues(lmat);
	//cout<<"Computed Fast EigenValues"<<endl;
	VectorC val = FastEigenVectors(lmat);
	//cout<<"Computed Eigenvalues and EigenVectors "<<endl;
	lmat.NormaliseRows();
	//cout<<"Computed Normalised Rows for EigenVectors"<<endl;
	MatrixC out = lmat.SubMatrix(lmat.Rows(), numclust);
	return out;
}

DListC<UIntT> SpectralClustering::GetKMeansClassification(const MatrixC &emat)
{
	//Convert the matrix into a set of SampleC<VectorC> objects
	SampleC<VectorC> sample;
	//Now create samples from each row of the eigenvector matrix
	for(UIntT i = 0; i < emat.Rows(); i ++)
	{
		VectorC vec(numclust); 
		vec.Fill(0.0);
		for(UIntT j = 0; j < numclust; j++)
		{
			vec[j] = emat[i][j];	
		}
		sample.Append(vec.Copy());
	}
	//Given a sample of data, now we can design a K-Means classifier
	DesignKMeansC km(numclust);
	ClassifierC kmcl = km.Apply(sample);
	//Given the classifier, now simply classify the data and output the result into a DListC<UIntT>
	DListC<UIntT> out;
	for(DArray1dIterC<VectorC> it(sample.DArray()); it; it++)
	{
		out.Append(kmcl.Classify((*it)));
	}
	return out;
}

ImageC<UIntT> SpectralClustering::Apply()
{
	//cout<<"Inside SpectralClustering Apply()"<<endl;
	Tuple2C<MatrixC,MatrixC> sd = ComputeSDMatrix(SHINCUT);
	//cout<<"Computed SD Matrix - "<<sd.Data1().Size()<<"\t"<<sd.Data2().Size()<<endl;
	MatrixC l = ComputeLaplacian(sd,UNNORMALISED);
	//cout<<"Computed Laplacian - "<<l.Size()<<endl;
	MatrixC e = GetNormalisedEigenVectors(l);
	//cout<<"Computed Eigenvectors - "<<e.Size()<<endl;
	DListC<UIntT> res = GetKMeansClassification(e);
	//cout<<"Got Image Classifications Using K-Means "<<res.Size()<<endl;
	ImageC<UIntT> diagimg(imrect,0);
	UIntT ind = 0;
	for(Array2dIterC<UIntT> it(diagimg); it; it++)
	{
		(*it) = res.Nth(ind);
		ind++;
	}
	return diagimg;
}
