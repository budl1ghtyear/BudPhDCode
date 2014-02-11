#include "LoadingFunctions.hh"
//Method if we want our original slow output
DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > MCDClustering(const FilenameC &im, const UIntT &c_type, const UIntT &im_type, const ImageRectangleC &irec)
{
	DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > result;
	switch(c_type)
	{
		case 0:
		{
			RealT h = 0.0;
			cout<<"You have chosen to use the Cascade MCD method, please enter the value of h you would like"<<endl;
			cin>>h;
			CascadedMCD mcd(LoadImage(im, im_type,irec),h);//load with image type we would like
			//Get the data we want to cluster and input into the clustering algorithm
			result = mcd.Apply();
			break;
		}
		case 1:
		{
			RealT h = 0.75;
			KSmirnovMCD mcd(LoadImage(im, im_type,irec),h); //load it up with the image type we would like
			//Get the data we want to cluster and input into the clustering algorithm
			result = mcd.Apply();
			break;	
		}
		default:
		{
			RealT h = 0.0;
			cout<<"You have chosen to use the Cascade MCD method, please enter the value of h you would like"<<endl;
			cin>>h;
			CascadedMCD mcd(LoadImage(im, im_type,irec),h);//load with image type we would like
			//Get the data we want to cluster and input into the clustering algorithm
			result = mcd.Apply();
			break;
		
		}
	}
	return result;
}

//Method for fast output
DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > MCDClustering(const ImageC<RealRGBValueC> &im, const UIntT &c_type, const UIntT &im_type, const ImageRectangleC &irec)
{
	DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > result;
	switch(c_type)
	{
		case 0:
		{
			RealT h = 0.0;
			cout<<"You have chosen to use the Cascade MCD method, please enter the value of h you would like"<<endl;
			cin>>h;
			CascadedMCD mcd(LoadImage(im, im_type,irec),h);//load with image type we would like
			//Get the data we want to cluster and input into the clustering algorithm
			result = mcd.Apply();
			break;
		}
		case 1:
		{
			RealT h = 0.75;
			KSmirnovMCD mcd(LoadImage(im, im_type,irec),h); //load it up with the image type we would like
			//Get the data we want to cluster and input into the clustering algorithm
			result = mcd.Apply();
			break;	
		}
		default:
		{
			RealT h = 0.0;
			cout<<"You have chosen to use the Cascade MCD method, please enter the value of h you would like"<<endl;
			cin>>h;
			CascadedMCD mcd(LoadImage(im, im_type,irec),h);//load with image type we would like
			//Get the data we want to cluster and input into the clustering algorithm
			result = mcd.Apply();
			break;
		
		}
	}
	return result;
}

