#include "PCAProjection.hh"

PCAProjectionC::PCAProjectionC(StringC const &pcadir, UIntT const &numcp)
{
	if(pcadir.Size() == 0)
	{
		cout<<"Directory was not instantiated"<<endl;
		directory="/vol/vssp/lip-tracking/StatisticalShapeModelsOfLips/PCA/";
	}
	else
		directory = const_cast<StringC &>(pcadir);
	
	FilenameC pmean = directory+"Mean.tsv";
	VectorC pca_m;
	FilenameC psdev = directory+"Eigenvalues.tsv";
	VectorC sdev;
	FilenameC ploadings = directory+"Eigenvectors.tsv";
	MatrixC pca_l;
	if((!LoadFileFromStream(pmean, pca_m))||(!LoadFileFromStream(psdev,sdev))||(!LoadFileFromStream(ploadings,pca_l)))
	{
		cerr<<"Could not load the PCA files"<<endl;
		exit(1);
	}
	//Call the remaining constructors
	pca_mean = pca_m.Copy();
	pca_loadings = pca_l.Copy();
	numcomp = numcp;
	pca_sdev = sdev.Copy();
	//cout<<"PCAProjection Constructor Called"<<endl;
}
PCAProjectionC::PCAProjectionC(VectorC const &mean, MatrixC const &eigenvect, VectorC const &eigenval, UIntT const &numcp):pca_mean(mean),pca_loadings(eigenvect),pca_sdev(eigenval)//IMPORTANT - THIS IMPLEMENTATION ASSUMES THAT THE EIGENVECTORS ARE IN EACH COLUMNS
{
	directory = "/vol/vssp/lip-tracking/StatisticalShapeModelsOfLips/PCA/";
	numcomp = numcp;
	//cout<<"PCAProjection Constructor Called"<<endl;
}
VectorC PCAProjectionC::PerformPCAProjection(SArray1dC<Point2dC> const &points)
{
	VectorC vect = PointsToVector(points);
	//Create the appropriate submatrix
	MatrixC submat = this->pca_loadings.SubMatrix(pca_loadings.Rows(),numcomp);
	FuncMeanProjectionC proj(pca_mean,submat.T()); //NOTE THAT this submat is transposed because the RAVL class expects a Row based eigenvector and we have a column based vector
	return proj.Apply(vect);
}
VectorC PCAProjectionC::PerformPCAProjection(Array1dC<Point2dC> const &points)
{
	return VectorC(PerformPCAProjection(ArrayToSArray_B(points)).Copy());
}
SArray1dC<Point2dC> PCAProjectionC::PerformInverseProjection(VectorC const &vect)
{
	MatrixC submat = this->pca_loadings.SubMatrix(pca_loadings.Rows(),numcomp);
	FuncLinearC proj(submat,this->pca_mean);
	VectorC vector = proj.Apply(vect);
	return SArray1dC<Point2dC>(VectorToPoints(vector).Copy());
}
