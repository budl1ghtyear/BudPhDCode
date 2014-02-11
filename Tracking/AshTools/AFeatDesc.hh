#ifndef VBAPAT_AFEATDESC_HEADER
#define VBAPAT_AFEATDESC_HEADER 1

#include "./AFeature2.hh"

//
// AFeature2C represents an affine-covariant feature with descriptor
//

namespace vbapAT {

   using namespace RavlN;

   class AFeatDescC : public AFeature2C
   {
   public:
      //:-
      // CONSTRUCTION/DESTRUCTION /////////////////////////////////////////////

      AFeatDescC(IntT length = 128);
      //: Default constructor

      AFeatDescC(RealT u,
                 RealT v,
                 VectorC desc);
      //: Constructor

      void Load(istream &is, Index2dC topLeft = Index2dC(0,0));
      //: Load from stream.

      void Save(ostream &os, Index2dC topLeft = Index2dC(0,0)) const;
      //: Save to stream.

      inline IntT Length() const
      {
         return m_descriptor.Size();
      }
      //: Return descriptor length

      inline const VectorC& Descriptor() const
      {
         return m_descriptor;
      }
      //: Return descriptor

      RealT EuclidDistance(const AFeatDescC& feature) const;
      //: Return Euclidian distance between feature descriptors

   protected:
      VectorC m_descriptor;
   };

}

#endif 
