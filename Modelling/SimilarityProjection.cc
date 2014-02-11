#include "SimilarityProjection.hh"

SimilarityProjectionC::SimilarityProjectionC(StringC const &dir)
{
	if(dir.Size() == 0)
	{
		cout<<"Directory was not instantiated"<<endl;
		directory = "/vol/vssp/lip-tracking/StatisticalShapeModelsOfLips/PCA/";
	}
	else
		directory = const_cast<StringC &>(dir);
	FilenameC affinemean = dir+"MeanPoints.txt";
	SArray1dC<Point2dC> aff_mean;
	if(!LoadFileFromStream(affinemean, aff_mean))
	{
		cerr<<"Could not load affine mean"<<endl;	
		exit(1);
	}
	a_mean = aff_mean.Copy();
	//cout<<"SimilarityProjection Constructor Called"<<endl;
}

SimilarityProjectionC::SimilarityProjectionC(SArray1dC<Point2dC> const &affinemean)
{
	directory = "/vol/vssp/lip-tracking/StatisticalShapeModelsOfLips/PCA/";
	a_mean = affinemean.Copy();
	//cout<<"SimilarityProjection Constructor Called"<<endl;
}

Tuple2C<SArray1dC<Point2dC>,VectorC> SimilarityProjectionC::ProjectToAffineSpace(SArray1dC<Point2dC> const &pts)
//This method will compute the similarity projection between pts and a_mean and project the pts points onto a_mean
//The returned points contain the residual of the two 
{
	 Affine2dC sim;
	 FitSimilarity(pts,a_mean,sim);
	 //Now project pts onto a_mean
	 //~ cout<<"Points inside similarity transformation: \n"<<pts<<endl;
	 SArray1dC<Point2dC> ptsinaffinespace = this->ProjectAffinePoints(pts,sim);
	 //~ cout<<"Points in affine space: \n"<<ptsinaffinespace<<endl;
	 //~ cout<<"Affine mean: \n"<<a_mean<<endl;
	 //Parameter composition - as written in AAMRotationScaleNormalisation
	 VectorC fixedParams(4);
	 Point2dC translation;Vector2dC scaling;RealT skew;RealT angle;
	 sim.Decompose(translation,scaling,skew,angle);
	 //~ cout<<scaling<<endl;
	 //~ cout<<"Applied affine transform is - "<<sim<<endl;
	 fixedParams[0] = translation.Row();
	 fixedParams[1] = translation.Col();
	 RealT scale = scaling[0];
	 fixedParams[2] = scale * cos(angle) - 1.0;
	 fixedParams[3] = scale * sin(angle);
	 //return Tuple2C<SArray1dC<Point2dC>,VectorC>(SubtractSArrays(ptsinaffinespace,a_mean),fixedParams);
	 return Tuple2C<SArray1dC<Point2dC>,VectorC>(ptsinaffinespace,fixedParams);
}
SArray1dC<Point2dC> SimilarityProjectionC::ProjectFromAffineSpace(VectorC const &vec, VectorC const &aff)
{
	SArray1dC<Point2dC> pts = VectorToPoints(vec);
	Vector2dC mean(aff[0],aff[1]);
    Matrix2dC sr(1.0+aff[2],aff[3],-aff[3],1.0+aff[2]);	
    Affine2dC sim(sr,mean);
    //~ cout<<"Input affine transform is : "<<sim<<endl;
	//~ cout<<"Applied affine transform is : "<<sim.Inverse()<<endl;
	//cout<<"Resulting vect in aff space - "<<PointsToVector(AddSArrays(pts,a_mean))<<endl;
	//SArray1dC<Point2dC> residual = ProjectAffinePoints(AddSArrays(pts,a_mean),sim.Inverse());
	SArray1dC<Point2dC> residual = ProjectAffinePoints(pts,sim.Inverse());
    return (residual);
}
SArray1dC<Point2dC> SimilarityProjectionC::ProjectAffinePoints(SArray1dC<Point2dC> const &pts, Affine2dC const &proj)
{
	SArray1dC<Point2dC> res(pts.Size());
	for(UIntT i = 0; i < pts.Size(); i++)
	{
		Vector2dC vec = proj * Vector2dC(pts[i].Row(),pts[i].Col());
		Point2dC pt(vec[0],vec[1]);
		res[i] = pt.Copy();
	}
	return res;
}
