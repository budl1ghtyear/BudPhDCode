#include "LipExtractorC.hh"

LipExtractorC::LipExtractorC(const ImageC<UIntT> &img, const ImageRectangleC &mR, const UIntT &num_pts)
{
	lip_img = img.Copy();
	num_points = num_pts;
	NUM_POINTS_EACH_SIDE = (num_points - 2) / 2;
	mouthRegion = mR;
	selected_Points = this->Apply(); //this will automatically instantiate the lip_bound_pts variable inside it
}

SArray1dC<Point2dC> LipExtractorC::Apply(void)
{
	ImageC<UIntT> boundary_img = PerformMorphology().Copy();
	//if(!Save("@X:Original Boundary Binary Lip Image",boundary_img)) cerr<<"Could not save binary lip image"<<endl;
	outerLipBoundary = BoundaryC(boundary_img,0);
	DListC<BoundaryC> outerLipBoundaryList = outerLipBoundary.OrderEdges();
	//cout<<"Number of boundary elements = "<<outerLipBoundaryList.Size()<<endl;
	UIntT maxLength=0;
	UIntT maxIndex=0;
	UIntT curIndex=0;
	for(DLIterC<BoundaryC> it(outerLipBoundaryList);it;it++)
	{
		//cout<<"Element is - "<<curIndex<<"\t with area = "<<Abs((*it).Area())<<endl;
		if (Abs((*it).Area()) > maxLength)
		{
			maxLength = Abs((*it).Area());
			maxIndex = curIndex;
			//cout<<"Max Index is Now - "<<maxIndex<<endl;
		}
		curIndex++;
		//cout<<"values are now - "<<curIndex<<"\t Max - "<<maxIndex<<endl;
	}
	outerLipBoundary = outerLipBoundaryList.Nth(maxIndex).Copy();
	//cout << "Selected outer lip boundary area: " << outerLipBoundary.Area() << endl;
	//cout << "Outer lip boundary pixels: " << outerLipBoundary << endl;
	lip_corner_points = GetLipCornerPoints();
	allVectors = GetLipNormalVectors();
	selected_Points = SelectLipBoundaryPoints().Copy();
	return selected_Points;
}

ImageC<UIntT> LipExtractorC::PerformMorphology(void)
{
	ImageC<UIntT> kernel(IndexRangeC(-1,1),IndexRangeC(-1,1));
	kernel.Fill(1);
	ImageC<UIntT> eroded_img,closed_img;
	MorphBinaryClose2d(lip_img,kernel,closed_img,(UIntT)1);
	BinaryErode(closed_img,kernel,eroded_img,(UIntT)1);
	ImageC<UIntT> boundary_img(lip_img.Frame(),0);
	for (Array2dIter3C<UIntT,UIntT,UIntT> it(boundary_img,closed_img,eroded_img,mouthRegion);it;it++)
	{
		it.Data1()=it.Data2()-it.Data3();
		if (it.Data1()>0)
		{
			SArray1dC<Point2dC> temp(1);
			temp[0]=Point2dC(it.Index().Row(),it.Index().Col()).Copy();
			lip_bound_pts.Append(temp);
		}
	}
	return boundary_img;
}

PairC<IndexC> LipExtractorC::GetLipCornerPoints(void)
{
	Point2dC centre;
	RealT major;
	RealT minor;
	RealT angle;
	FitEllipse(lip_bound_pts,lip_ellipse);
	lip_ellipse.EllipseParameters(centre, major, minor, angle);
	//cout << "Centre coords: " << centre << ", Major axis length: " << major << ", Minor axis length: " << minor << ", rotation angle: " << angle << endl;
	SArray1dC<RealT> major_axis_dist(lip_bound_pts.Size());
	major_axis_dist.Fill((RealT)100.0);
	SArray1dC<RealT> minor_axis_dist(lip_bound_pts.Size());
	minor_axis_dist.Fill((RealT)0.0);
	LineABC2dC major_line(-Sin(angle), Cos(angle), -(-Sin(angle)*centre.Row()+Cos(angle)*centre.Col()));
	LineABC2dC minor_line(Cos(angle), Sin(angle), -(Cos(angle)*centre.Row()+Sin(angle)*centre.Col()));
	for (SArray1dIter3C<Point2dC, RealT, RealT> it(lip_bound_pts, major_axis_dist, minor_axis_dist);it;it++)
	{
		it.Data2()=major_line.SqrEuclidDistance(it.Data1());
		it.Data3()=minor_line.SignedDistance(it.Data1());
	}
	IndexC leftBound = minor_axis_dist.IndexOfMin();
	//cout << "Leftmost point is " << lip_bound_pts[leftBound] << "at a distance of " << minor_axis_dist[leftBound] << endl;
	//cout << "Distance from major axis: " << major_axis_dist[leftBound]<<endl;
	if (major_axis_dist[leftBound] < 100)
	{
		//Bingo
		//cout << "Nice left side at " << lip_bound_pts[leftBound] << endl;
	}
	IndexC rightBound = minor_axis_dist.IndexOfMax();
	//cout << "Rightmost point is " << lip_bound_pts[rightBound] << "at a distance of " << minor_axis_dist[rightBound] << endl;
	//cout << "Distance from major axis: " << major_axis_dist[rightBound]<<endl;
	if (major_axis_dist[rightBound] < 100)
	{
		//Bingo
		//cout << "Nice right side at " << lip_bound_pts[rightBound] << endl;
	}
	return (PairC<IndexC> (leftBound, rightBound));	
}

SArray1dC<Vector2dC> LipExtractorC::GetLipNormalVectors(void)
{
	Point2dC centre;
	RealT major;
	RealT minor;
	RealT angle;
	lip_ellipse.EllipseParameters(centre, major, minor, angle);
	IndexC leftBound = lip_corner_points.A();
	IndexC rightBound = lip_corner_points.B();
	Vector2dC toLeft(lip_bound_pts[leftBound].Row()-centre.Row(), lip_bound_pts[leftBound].Col()-centre.Col());
	Vector2dC toRight(lip_bound_pts[rightBound].Row()-centre.Row(), lip_bound_pts[rightBound].Col()-centre.Col());
	
	RealT ellipseAngleCos = toRight.Dot(toLeft)/(toLeft.Norm()*toRight.Norm());
	//RealT ellipseAngleCos = toLeft.Dot(toRight)/(toLeft.Norm()*toRight.Norm());
	//cout << "Cosine for extremity vecor angle (top side): " << ellipseAngleCos << endl;
	RealT ellipseAngleSin = toRight.Cross(toLeft)/(toLeft.Norm()*toRight.Norm());
	//RealT ellipseAngleSin = toLeft.Cross(toRight)/(toLeft.Norm()*toRight.Norm());
	//cout << "Sine for extremity vecor angle (top side): " << ellipseAngleSin << endl;
	RealT angleTopSide = ACos(ellipseAngleCos);
	if (ellipseAngleSin<0) angleTopSide = (RealT)2.*pi - angleTopSide;
	//cout << "Vecor angle (top side): " << angleTopSide << endl;
	RealT topPartitioningAngle = angleTopSide/((RealT)NUM_POINTS_EACH_SIDE+(RealT)1.);
	//cout << "Interval angle (top side): " << topPartitioningAngle << endl;
	ComplexC leftEnd(toLeft.X(), toLeft.Y());
	ComplexC partitionMultiplierTop(cos(topPartitioningAngle), -sin(topPartitioningAngle));
	SArray1dC<ComplexC> orientationsTop(NUM_POINTS_EACH_SIDE+1);
	orientationsTop[0]=leftEnd;
	for (UIntT i=1;i<=NUM_POINTS_EACH_SIDE;i++)
	{
		orientationsTop[i]=orientationsTop[i-1]*partitionMultiplierTop;
	}
	//cout << "Top vectors:" << endl << orientationsTop << endl;

	RealT bottomPartitioningAngle = ((RealT)2.*pi-angleTopSide)/((RealT)NUM_POINTS_EACH_SIDE + (RealT)1.);
	//cout << "Interval angle (bottom side): " << bottomPartitioningAngle << endl;
	ComplexC rightEnd(toRight.X(), toRight.Y());
	ComplexC partitionMultiplierBottom(cos(bottomPartitioningAngle), -sin(bottomPartitioningAngle));
	SArray1dC<ComplexC> orientationsBottom(NUM_POINTS_EACH_SIDE+1);
	orientationsBottom[0]=rightEnd;
	for (UIntT i=1;i<=NUM_POINTS_EACH_SIDE;i++)
	{
		orientationsBottom[i]=orientationsBottom[i-1]*partitionMultiplierBottom;
	}
	//cout << "Bottom vectors:" << endl << orientationsBottom << endl;
	//~ SArray1dC<ComplexC> loopCloser(1);
	//~ loopCloser[0]=leftEnd;
	SArray1dC<ComplexC> allOrientations;
	allOrientations.Append(orientationsTop);
	allOrientations.Append(orientationsBottom);
	//allOrientations.Append(loopCloser); //THIS LINE CAUSES THE FIRST POINT TO BE REPEATED. IDEALLY WE WANT A CLOSED BSPLINE IMPLEMENTATION BUT TILL THEN, THIS WILL HAVE TO DO!!
	cout << "Orientations: " << endl << allOrientations <<endl;
	SArray1dC<Vector2dC> allVectors(allOrientations.Size());
	for (SArray1dIter2C<Vector2dC,ComplexC> it(allVectors,allOrientations);it;it++)
	{
		it.Data1()=Vector2dC(it.Data2().Re(), it.Data2().Im());
		//DrawLine(boundary_img,(UIntT)127,(UIntT)127,Index2dC(centre.Row(),centre.Col()),Index2dC(centre.Row()+it.Data1().X(),centre.Col()+it.Data1().Y()));
	}
	return allVectors;
}

SArray1dC<Point2dC> LipExtractorC::SelectLipBoundaryPoints(void)
{
	Point2dC centre;
	RealT major;
	RealT minor;
	RealT angle;
	lip_ellipse.EllipseParameters(centre, major, minor, angle);
	SArray1dC<LineABC2dC> normalLines(allVectors.Size());
	SArray1dC<Point2dC> selectedPoints(allVectors.Size());
	for (SArray1dIter3C<Point2dC,LineABC2dC,Vector2dC> it(selectedPoints,normalLines,allVectors);it;it++)
	{
		RealT sign;
		if ((it.Index()>0) && (it.Index()<=NUM_POINTS_EACH_SIDE))
		{
			sign=1;
		}
		else if ((it.Index()>NUM_POINTS_EACH_SIDE+1) && (it.Index()<=2*NUM_POINTS_EACH_SIDE+1))
		{
			sign=-1;
		}
		else
		{	
			sign=0;
		}
		//it.Data2()=LineABC2dC(-it.Data3().Y(), -it.Data3().X(), centre.Row()*(centre.Col()+it.Data3().Y())-(centre.Row()+it.Data3().X())*centre.Col());
		it.Data2()=LineABC2dC(-it.Data3().Y(), -it.Data3().X(), centre.Row()*it.Data3().Y()+it.Data3().X()*centre.Col());
		//cout << "Line equation : " << it.Data2().A() << "*x + " << it.Data2().B() << "*y + " << it.Data2().C() << " = 0" << endl;
		RealT minDist = (RealT)100.;
		IndexC minIndex = 0;
		for (SArray1dIterC<Point2dC> it2(lip_bound_pts);it2;it2++)
		{
			bool pointAllowed = false;
			for (DLIterC<CrackC> it3(outerLipBoundary);it3;it3++)
			{
				if ((*it3).RPixel()==Index2dC((*it2).X(), (*it2).Y()))
				{
					pointAllowed = true;
					break;
				}
			}
			if (pointAllowed)
			{	if (it.Data2().Distance(*it2) < minDist)
				{
					if (sign == 0)
					{
						if (it.Index()==NUM_POINTS_EACH_SIDE+1)
						{
							if ((*it2).Y()-centre.Col()>0)
							{
								minDist = it.Data2().Distance(*it2);
								minIndex = it2.Index();
								//cout << "New minimum of " << minDist << " at " << minIndex << endl;
							}
						}
						else //(it.Index()==0) || (it.Index()==2*NUM_POINTS_EACH_SIDE+1)
						{
							if ((*it2).Y()-centre.Col()<0)
							{
								minDist = it.Data2().Distance(*it2);
								minIndex = it2.Index();
								//cout << "New minimum of " << minDist << " at " << minIndex << endl;
							}
						}
					}
					else
					{
						if (sign*((*it2).X()-centre.Row())>=0)
						{
							minDist = it.Data2().Distance(*it2);
							minIndex = it2.Index();
							//cout << "New minimum of " << minDist << " at " << minIndex << endl;
						}
					}
				}
			}
		}
		it.Data1()=lip_bound_pts[minIndex];
	}
	for (SArray1dIterC<Point2dC> it(selectedPoints);it;it++)
	{
		//boundary_img
		DrawCross(lip_img,(UIntT)127,Index2dC((*it).Row(),(*it).Col()),3);
	}
	cout << "Contour points to be used: " << endl << selectedPoints <<endl;
	if(!Save("@X:Boundary Binary Lip Image",lip_img)) cerr<<"Could not save binary lip image"<<endl;
	return selectedPoints;
}
