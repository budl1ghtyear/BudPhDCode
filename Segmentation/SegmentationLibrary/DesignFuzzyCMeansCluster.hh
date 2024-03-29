// This file is part of RAVL, Recognition And Vision Library 
// Copyright (C) 2003, University of Surrey
// This code may be redistributed under the terms of the GNU Lesser
// General Public License (LGPL). See the lgpl.licence file for details or
// see http://www.gnu.org/copyleft/lesser.html
// file-header-ends-here
#ifndef RAVL_MYDESIGNFUZZYCMEANSCLUSTER_HEADER
#define RAVL_MYDESIGNFUZZYCMEANSCLUSTER_HEADER
////////////////////////////////////////////////////////////////////////////
//! author="Robert Crida"
//! lib=RavlPatternRec
//! date="30/1/1998"
//! userlevel=Develop
//! file="Ravl/PatternRec/Cluster/DesignFuzzyCMeansCluster.hh"
//! docentry="Ravl.API.Pattern Recognition.Cluster;Ravl.API.Pattern Recognition.Classifier.DesignClassifier"
//! rcsid="$Id: DesignFuzzyCMeansCluster.hh 5240 2005-12-06 17:16:50Z plugger $"
//! example="exFuzzyCMeansCluster.cc"

#include "Ravl/PatternRec/DesignCluster.hh"
#include "Ravl/PatternRec/DistanceSqrEuclidean.hh"

#include "Ravl/PatternRec/ClassifierFuzzyCMeans.hh"
#include "Ravl/StdMath.hh"
#include "Ravl/Random.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/PatternRec/SampleIter.hh"

using namespace RavlN;

//namespace RavlN {

  //! userlevel=Develop
  //: Fuzzy C-means clustering implementation class.
  //
  // Implementation class for fuzzy C-means clustering algorithm. The
  // DesignFuzzyCMeansClusterC handle class should be used.
  
  class MyDesignFuzzyCMeansClusterBodyC : public DesignClusterBodyC
  {
  public:
    MyDesignFuzzyCMeansClusterBodyC (UIntT numClasses, const DistanceC &distance, RealT m);
    //: Constructor
    //!param: numClasses - specifies how many clusters to assume
    //!param: distance   - metric to use when measuring sample distance from class centres
    //!param: m          - class membership weighting exponent

    MyDesignFuzzyCMeansClusterBodyC (istream &in);
    //: Constructs from stream

    MyDesignFuzzyCMeansClusterBodyC (const MyDesignFuzzyCMeansClusterBodyC &oth);
    //: Copy Constructor

    virtual RCBodyVC & Copy () const;
    //: Makes a deep copy and is virtual
    
    virtual void Initialise (UIntT numClasses, const SampleC<VectorC> &trainX);
    //: Initialises the clusters
    //!param: numClasses - number of clusters that will be set up
    //!param: trainX     - training samples used to determine initial centers
    
    virtual FunctionC Apply(const SampleC<VectorC> &trainX);
    //: Used to learn relationships c=f(X) and Y=f(X)
    
    virtual bool Save (ostream &out) const;
    //: Writes object to stream, can be loaded using constructor

    friend class DesignFuzzyCMeansClusterC;
    //: Handle class
    
  protected:
    UIntT sizeY;
    UIntT sizeX;
    
    SArray1dC<VectorC> _centres;
    //: Cluster centres

    DistanceC _distance;
    //: Metric for measuring distance between samples and cluster centres
    
    RealT _m;
    //: Class membership weighting exponent

    bool _initialised;
  };

  //! userlevel=Normal
  //: Fuzzy C-means clustering algorithm
  //
  // Unsupervised classification algorithm which performs fuzzy C-means
  // clustering. It is necessary to specify the number of classes that are
  // expected. Initial class centres are determined by assigning each sample
  // random membership to the class centres.

  class MyDesignFuzzyCMeansClusterC
    : public DesignClusterC
  {
  public:
    MyDesignFuzzyCMeansClusterC (UIntT numClasses, const DistanceC &distance = DistanceSqrEuclideanC(), RealT m=2.0)
      : DesignClusterC(*(new MyDesignFuzzyCMeansClusterBodyC (numClasses,distance,m))) 
    {}
    //: Constructor
    //!param: numClasses - number of classes or clusters to look for
    //!param: distance   - metric to use when measuring sample distance from class centres
    //!param: m          - class membership weighting exponent
  };
//}

#endif
