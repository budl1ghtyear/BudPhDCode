#include "ModellingUtilityFunctions.hh"

DListC<SArray1dC<Point2dC> > ModellingUtilityFunctions::LoadTRESData(const FilenameC& f)
//: Generate a data structure that contains the tres specified lip-shape co-ordinates
{
	//LOAD THE LANDMARKED DATA
	DListC<SArray1dC<Point2dC> > lp_coords;
	IStreamC is(f);
	while(!is.IsEndOfStream())
	{
		SArray1dC<Point2dC> pts;
		is >> pts;
		if(pts.Size() > 0)
		{
			//If data used is EJ, need to reverse order
			SArray1dC<Point2dC> to_append = ReorderPoints(pts);
			lp_coords.Append(to_append);			
		}				
	}	
	return lp_coords;
}

SArray1dC<Point2dC> ModellingUtilityFunctions::ReorderPoints(const SArray1dC<Point2dC> &lp_pts)
//:This method is specifically applicable to the data provided by EJ using LP to landmark lip points
//:Since his method uses the CImg library, his x and y coordinates are reversed
//:EJ data is loaded straight into a Point2dC object, therefore, Old.Row() = New.Col() and vice versa.
{
	SArray1dC<Point2dC> new_lp(lp_pts.Size());
	for(UIntT i = 0; i < lp_pts.Size(); i++)
	{
		new_lp[i] = Point2dC(lp_pts[i].Col(), lp_pts[i].Row());
	}
	return new_lp;
}
 
Array1dC<Point2dC> ModellingUtilityFunctions::SArrayToArray(const SArray1dC<Point2dC> &pts)
//:Data structure converter
{
	Array1dC<Point2dC> res(pts.Size());
	for(UIntT i = 0; i < pts.Size(); i++)
	{
		res[i] = pts[i].Copy();
	}
	return res;
}

ImageC<RealRGBValueC> ModellingUtilityFunctions::LoadImage(const FilenameC& f)
//:Load Image
{
	ImageC<RealRGBValueC> img;
	if(!Load(f,img)) cerr<<"Could not load the image file "<<f<<endl;
	return img;
}

VectorC ModellingUtilityFunctions::ArrayToVector(const Array1dC<Point2dC> &arr)
//:Convert the control point array to a vector - {X0....Xm,Y0...Ym} i.e. {Col[0]....Col[m],Row[0]...Row[m]}
{
	VectorC vec(arr.Size()*2);
	vec.Fill(0);
	for(UIntT i = 0 ; i < arr.Size() ; i++)
	{
		vec[i] = arr[i].Col();
		vec[i + arr.Size()] = arr[i].Row();
	}
	return vec;
}

Point2dC ModellingUtilityFunctions::TranslationNormalise(const Array1dC<Point2dC> &arr)
//:Perform translation normalisation
{
	Point2dC mean = ComputeMeanPoint(arr);	
	//Now subtract the mean position from each element inside the Array1dC
	//Array1dC<Point2dC> res = arr - mean;
	return mean;
}

RealT ModellingUtilityFunctions::RotationNormalise(const Array1dC<Point2dC> &arr)
//:	Perform rotation and translation normalisation using the lip corners
{
	//Lip Corners 
	Tuple2C<Point2dC,Point2dC> lip_corners(arr[0].Copy(),arr[8].Copy());
	//Horizontal x-axis
	Point2dC lpt(lip_corners.Data1().Copy());Point2dC rpt(lip_corners.Data1().Row(),lip_corners.Data2().Col());
	//Compute rotation angle
	Vector2dC x_axis = rpt - lpt;
	Vector2dC lip_crn = lip_corners.Data2() - lip_corners.Data1();
	RealT rot = lip_crn.Angle() - x_axis.Angle();    
	return rot;
}

Vector2dC ModellingUtilityFunctions::ScaleNormalise(const RealT &mn, const RealT &md)
{
	Vector2dC res((RealT)(mn/md),(RealT)(mn/md));
	return res;
}
Array1dC<Point2dC> ModellingUtilityFunctions::PerformProjection(const Affine2dC &aff, const Array1dC<Point2dC> &pts)
{
	Array1dC<Point2dC> res(pts.Size());
	for(UIntT i = 0; i < pts.Size(); i++)
	{
		Vector2dC vec(pts[i].Row(),pts[i].Col());
		Vector2dC proj_vec = aff * vec;
		Point2dC proj_pt(proj_vec.Row(),proj_vec.Col());
		res[i] = proj_pt.Copy();
	}
	return res;
}

Array1dC<Point2dC> ModellingUtilityFunctions::AffineNormalise(const Array1dC<Point2dC> &pts, const Point2dC &trans, const Vector2dC &scale,const RealT &rotation)
{
	Affine2dC aff(scale,rotation,trans);
	//cout<<"Input parameters - "<<trans<<"\t scale "<<scale<<"\t rotation - "<<rotation<<endl;
	//cout<<"Input array to Affine Normalise "<<pts<<endl;
	//cout<<"Affine transform "<<aff.Translation()<<"\n SR Matrix "<<aff.SRMatrix()<<endl;
	Array1dC<Point2dC> res = PerformProjection(aff,pts);
	//cout<<"Output array to Affine Normalise "<<res<<endl;
	return res;
}

Point2dC ModellingUtilityFunctions::ComputeMeanPoint(const Array1dC<Point2dC> &pts)
{
		Point2dC mean = pts.Sum();
		mean = mean / (RealT)(pts.Size());
		return mean;
}
SArray1dC<Point2dC> ModellingUtilityFunctions::ArrayToSArray(const Array1dC<Point2dC> &pts)
{
	SArray1dC<Point2dC> res(pts.Size());
	for(UIntT i = 0; i < pts.Size(); i++)
	{
		res[i] = pts[i].Copy();
	}
	return res;	
}
