// Created on 21/01/2009.
// Modified on 21/01/2009.
// Ravl includes

#include "Ravl/IO.hh"
#include "Ravl/DP/FileFormatIO.hh"
#include "Ravl/Array1dIter.hh"

#include "./AFeatDesc.hh"

namespace vbapAT
{

   //: Default constructor
   AFeatDescC::AFeatDescC(IntT length)
      : AFeature2C(),
        m_descriptor(length)
   {
   }

   //: Constructor
   AFeatDescC::AFeatDescC(RealT u,
                          RealT v,
                          VectorC desc)
      : AFeature2C(u, v),
        m_descriptor(desc)
   {
   }

   //: Return Euclidian distance between feature descriptors
   RealT AFeatDescC::EuclidDistance(const AFeatDescC& feature) const
   {
      return m_descriptor.EuclidDistance(feature.Descriptor());
   }

   void AFeatDescC::Load(istream &is, Index2dC topLeft)
   {
      is >> m_u >> m_v;
      m_u += topLeft.Col();
      m_v += topLeft.Row();
      for(Array1dIterC<RealT> it(m_descriptor); it; it++) {
        is >> *it;
      }
   }

   void AFeatDescC::Save(ostream &os, Index2dC topLeft) const
   {
      os << m_u-topLeft.Col() <<  " " << m_v-topLeft.Row();
      for(Array1dIterC<RealT> it(m_descriptor); it; it++) {
        os << " " << *it;
      }
   }

}
