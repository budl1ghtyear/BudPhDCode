
template <class T>
LipClusteringC<T>::LipClusteringC(const ImageC<T> &img, const RealT &hval = (RealT)0.75, const UIntT &numclust = (UIntT) 3):h(hval),n(numclust)
{
	inp_img = img.Copy();
}

template <class T>
DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > LipClusteringC<T>::Apply(const ClusteringType &ctype)
{
	DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > result;
	switch(ctype)
	{
		case CMCD:
		{
			CascadedMCD mcd(inp_img,h);
			//Get the data we want to cluster and input into the clustering algorithm
			result = mcd.Apply();
			break;
		}
		case KSMCD:
		{
			KSmirnovMCD mcd(inp_img,h); //load it up with the image type we would like
			//Get the data we want to cluster and input into the clustering algorithm
			result = mcd.Apply();
			break;	
		}
		case KM:
		{	
			DesignKMeansC km(n);
			SampleC<VectorC> smpl = GetSample(inp_img);
			ClassifierC cl = km.Apply(smpl);
			ImageC<UIntT> clustimg(inp_img.Frame(),0);
			UIntT dim = inp_img[0][0].Size();
			DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > out;
			SArray1dC<DListC<Tuple2C<VectorC,Index2dC> > > pix(n);
			for(Array2dIterC<T> it(inp_img); it; it++)
			{
				VectorC v(dim);
				for(UIntT i = 0; i < dim; i++)
				{
					v[i] = (*it)[i];
				}
				UIntT label = cl.Classify(v);
				clustimg[it.Index()] = label;
				//Insert VectorC and Index into the output cluster
				Tuple2C<VectorC,Index2dC> tpl(v.Copy(), it.Index());
				pix[label].Append(tpl);
			}
			//Now obtain the MeanCovariance Values for each cluster
			for(SArray1dIterC<DListC<Tuple2C<VectorC,Index2dC> > > it(pix); it; it++)
			{
				SumsNd2C sum(dim);
				bool SampleStats = false;
				for(DLIterC<Tuple2C<VectorC,Index2dC> > it2((*it)); it2; it2++)
				{
					sum += (*it2).Data1();
				}
				Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > elem(sum.MeanCovariance(SampleStats).Copy(),(*it));
				out.Append(elem.Copy());
			}
			result = out.Copy();	
			break;		
		}
		case FCM:
		{
			DesignFuzzyCMeansClusterC km(n);
			SampleC<VectorC> smpl = GetSample(inp_img);
			ClassifierC cl = km.Apply(smpl);
			ImageC<UIntT> clustimg(inp_img.Frame(),0);
			UIntT dim = inp_img[0][0].Size();
			DListC<Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > > out;
			SArray1dC<DListC<Tuple2C<VectorC,Index2dC> > > pix(n);
			for(Array2dIterC<T> it(inp_img); it; it++)
			{
				VectorC v(dim);
				for(UIntT i = 0; i < dim; i++)
				{
					v[i] = (*it)[i];
				}
				UIntT label = cl.Classify(v);
				clustimg[it.Index()] = label;
				//Insert VectorC and Index into the output cluster
				Tuple2C<VectorC,Index2dC> tpl(v.Copy(), it.Index());
				pix[label].Append(tpl);
			}
			//Now obtain the MeanCovariance Values for each cluster
			for(SArray1dIterC<DListC<Tuple2C<VectorC,Index2dC> > > it(pix); it; it++)
			{
				SumsNd2C sum(dim);
				bool SampleStats = false;
				for(DLIterC<Tuple2C<VectorC,Index2dC> > it2((*it)); it2; it2++)
				{
					sum += (*it2).Data1();
				}
				Tuple2C<MeanCovarianceC,DListC<Tuple2C<VectorC,Index2dC> > > elem(sum.MeanCovariance(SampleStats).Copy(),(*it));
				out.Append(elem.Copy());
			}
			result = out.Copy();
			break;					
		}
		default:
		{
			result = Apply(CMCD) ;
			break;
		}
	}
	return result;	
}
