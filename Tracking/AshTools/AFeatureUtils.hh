// Modified 23/01/2009.

#ifndef VBAPAT_AFEATUREUTILS_HEADER
#define VBAPAT_AFEATUREUTILS_HEADER 1

#include "Ravl/Stream.hh"
#include "Ravl/Vector.hh"
#include "Ravl/Ellipse2d.hh"
#include "Ravl/IndexRange2d.hh"
#include "Ravl/String.hh"
#include "Ravl/DList.hh"
// Added
#include "Ravl/DLIter.hh"
// Added

#include "./AFeature2.hh"
#include "./AFeatDesc.hh"

//
// AFeature2C represents an affine-covariant feature without descriptor
//

namespace vbapAT {

   using namespace RavlN;

   void Load(const StringC& filename, DListC<AFeature2C>& lp, const Index2dC& topLeft)
   {
      IStreamC is(filename);
      RealT version;
      is >> version;
      IntT nb_features;
      is >> nb_features;
      for(IntT featInd=0; featInd<nb_features; featInd++) {
         AFeature2C feature;
         feature.Load(is, topLeft);
         lp+=feature;
      }
   }
   // Reads a list of feature points from a file

   void Load(const StringC& filename, DListC<AFeatDescC>& lp, const Index2dC& topLeft)
   {
      IStreamC is(filename);
      IntT length;
      is >> length;
      IntT nb_features;
      is >> nb_features;
      for(IntT featInd=0; featInd<nb_features; featInd++) {
         AFeatDescC feature;
         feature.Load(is, topLeft);
         lp+=feature;
      }
   }
   // Reads a list of feature points (with descriptors) from a file

/*
   void Load2(const StringC& filename, DListC<AFeature2C>& lp, const Index2dC& topLeft)
   {
      IStreamC is(filename);
      IntT length;
      is >> length;
      IntT nb_features;
      is >> nb_features;
      for(IntT featInd=0; featInd<nb_features; featInd++) {
         AFeature2C feature;
         feature.Load2(is, topLeft);
         lp+=feature;
      }
   }
   // Reads a list of feature points (with descriptors) from a file
*/
/*   void Load(const StringC& filename, DListC<ACFeatureCorrespondenceC>& lp)
   {
      IStreamC is(filename);
      IntT nb_features;
      is >> nb_features;
      for(IntT featInd=0; featInd<nb_features; featInd++) {
         ACFeatureCorrespondenceC feature;
         is >> feature;
         lp+=feature;
      }
   }
   // Reads a list of feature correspondences from a file*/

   void Save(const StringC& filename, DListC<AFeature2C>& lp, const Index2dC& topLeft)
   {
      OStreamC os(filename);
      os << 1.0 << endl;
      os << lp.Size() << endl;
      for(DLIterC<AFeature2C> lfit(lp); lfit; lfit++) {
         lfit->Save(os, topLeft);
         os << "\n";
      }
   }
   // Save a list of feature points to a file

   void Save(const StringC& filename, DListC<AFeatDescC>& lp, const Index2dC& topLeft)
   {
      OStreamC os(filename);
      os << lp.First().Descriptor().Size() << endl;
      os << lp.Size() << endl;
      for(DLIterC<AFeatDescC> lfit(lp); lfit; lfit++) {
         lfit->Save(os, topLeft);
         os << "\n";
      }
   }
   // Save a list of feature points (with descriptors) to a file

/*   void Save(const StringC& filename, DListC<ACFeatureCorrespondenceC>& lp)
   {
      OStreamC os(filename);
      os << lp.Size() << endl;
      for(DLIterC<ACFeatureCorrespondenceC> lfit(lp); lfit; lfit++) {
         os << *lfit << endl;
      }
   }
   // Save a list of feature correspondences to a file*/

}

#endif 
