// Ravl includes

#include "Ravl/IO.hh"
#include "Ravl/DP/FileFormatIO.hh"

#include "./AFeature.hh"

namespace vbapAT
{

   AFeatureC::AFeatureC()
      : m_u(0.0),
        m_v(0.0),
        m_a(0.0),
        m_b(0.0),
        m_c(0.0)
   {
   }
   //: Default constructor

   AFeatureC::AFeatureC(RealT u,
                        RealT v,
                        RealT a,
                        RealT b,
                        RealT c)
      : m_u(u),
        m_v(v),
        m_a(a),
        m_b(b),
        m_c(c)
   {
   }
   //: Constructor

   //: Return Ellipse2dC representation of the feature
   Ellipse2dC AFeatureC::Ellipse() const
   {
      TFVectorC<RealT,6> coeffs(6);
      coeffs[0] = m_c;
      coeffs[1] = 2.0*m_b;
      coeffs[2] = m_a;
      coeffs[3] = -2.0*(m_b*m_u+m_c*m_v);
      coeffs[4] = -2.0*(m_a*m_u+m_b*m_v);
      coeffs[5] = m_a*Sqr(m_u)+m_c*Sqr(m_v)+2.0*m_b*m_u*m_v-1.0;
      return Ellipse2dC(coeffs);
   }

   //: Return a bounding box containing the feature
   IndexRange2dC AFeatureC::BoundingRectangle() const
   {
      return IndexRange2dC(Index2dC(m_v, m_u), Ceil(Height()), Ceil(Width()));
   }

   //: Is this pixel part of the feature?
   bool AFeatureC::Contains(const Index2dC& pix) const
   {
      RealT x = pix.Col()-m_u;
      RealT y = pix.Row()-m_v;
      return m_a*Sqr(x)+2.0*m_b*x*y+m_c*Sqr(y)<=1.0;
   }

   //: Is this image point part of the feature?
   bool AFeatureC::Contains(const Point2dC& pt) const
   {
      RealT x = pt.Col()-m_u;
      RealT y = pt.Row()-m_v;
      return m_a*Sqr(x)+2.0*m_b*x*y+m_c*Sqr(y)<=1.0;
   }

   void AFeatureC::Load(istream &is, Index2dC topLeft)
   {
      is >> m_u >> m_v >> m_a >> m_b >> m_c;
      m_u += topLeft.Col();
      m_v += topLeft.Row();
   }

   void AFeatureC::Save(ostream &os, Index2dC topLeft) const
   {
      os << m_u-topLeft.Col() <<  " " << m_v-topLeft.Row() << " " << m_a << " " << m_b << " " << m_c;
   }

}
