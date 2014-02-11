#ifndef VBAPAT_AFEATURE_HEADER
#define VBAPAT_AFEATURE_HEADER 1

#include "Ravl/Stream.hh"
#include "Ravl/Vector.hh"
#include "Ravl/Ellipse2d.hh"
#include "Ravl/IndexRange2d.hh"
#include "Ravl/DList.hh"

//
// AFeatureC represents an affine-covariant feature without descriptor
//

namespace vbapAT {

   using namespace RavlN;

   class AFeatureC
   {
   public:
      //:-
      // CONSTRUCTION/DESTRUCTION /////////////////////////////////////////////

      AFeatureC();
      //: Default constructor

      AFeatureC(RealT u,
                RealT v,
                RealT a,
                RealT b,
                RealT c);
      //: Constructor

      Ellipse2dC Ellipse() const;
      //: Return Ellipse2dC representation of the feature

      IndexRange2dC BoundingRectangle() const;
      //: Return a bounding box containing the feature

      bool Contains(const Index2dC& pix) const;
      //: Is this pixel part of the feature?

      bool Contains(const Point2dC& pt) const;
      //: Is this image point part of the feature?

      inline Point2dC Center() const
      {
         return Point2dC(m_v, m_u);
      }
      //: Return feature centre

      inline RealT Row() const
      {
        return m_v;
      }
      //: Return feature row

      inline RealT Col() const
      {
         return m_u;
      }
      //: Return feature column

      inline RealT Height() const
      {
         return 2*Sqrt(m_a/(m_a*m_c-Sqr(m_b)));
      }
      //: Return feature height

      inline RealT Width() const
      {
         return 2*Sqrt(m_c/(m_a*m_c-Sqr(m_b)));
      }
      //: Return feature width

      inline RealT MajorAxisLength() const
      {
         return 2*Sqrt(2/(m_a+m_c-Abs(m_a-m_c)*Sqrt(1+4*Sqr(m_b)/Sqr(m_a-m_c))));
      }
      //: Return major axis length

      inline RealT MinorAxisLength() const
      {
         return 2*Sqrt(2/(m_a+m_c+Abs(m_a-m_c)*Sqrt(1+4*Sqr(m_b)/Sqr(m_a-m_c))));
      }
      //: Return minor axis length

      inline RealT VerticalShift(const AFeatureC& feat) const
      {
         return m_v-feat.Row();
      }
      //: Return vertical shift between features

      void Load(istream &is, Index2dC topLeft = Index2dC(0,0));
      //: Load from stream.

      void Save(ostream &os, Index2dC topLeft = Index2dC(0,0)) const;
      //: Save to stream.

   protected:
      RealT m_u;
      RealT m_v;
      RealT m_a;
      RealT m_b;
      RealT m_c;
   };

}

#endif 
