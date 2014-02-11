#ifndef VBAPAT_AFEATURE_2_HEADER
#define VBAPAT_AFEATURE_2_HEADER 1

#include "Ravl/Stream.hh"
#include "Ravl/Vector.hh"
#include "Ravl/Ellipse2d.hh"
#include "Ravl/IndexRange2d.hh"
#include "Ravl/DList.hh"

#include "Ravl/SArray1d.hh"
#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ByteRGBValue.hh"
#include "Ravl/3D/PinholeCameraArray.hh"

#include "../MultiViewToolkit/Texture/MeshVisibility.hh"

//
// AFeatureC represents an affine-covariant feature without descriptor
// Modified on 06/01/2009 by Ashish Doshi.

namespace vbapAT {

   using namespace RavlN;
   using namespace Ravl3DN;
   using namespace RavlImageN;
   using namespace MVTTextureN;

   class AFeature2C
   {
   public:
      //:-
      // CONSTRUCTION/DESTRUCTION /////////////////////////////////////////////

      AFeature2C();
      //: Default constructor

      AFeature2C(RealT u,
                 RealT v);
      //: Constructor

//      AFeature2C::AFeature2C(RealT u1,
//                             RealT v1,
//			     RealT u2,
//			     RealT v2);
      //: Constructor

      AFeature2C(const TriMeshC& mesh, const PinholeCameraArrayC& cameras)
         : m_mesh(mesh), m_cameras(cameras)
      {}
      //: Construct from the mesh and camera calibration

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
/*
      inline RealT Row1() const
      {
         return m_v1;
      }
      //: Return feature row

      inline RealT Col1() const
      {
         return m_u1;
      }
      //: Return feature column

      inline RealT Row2() const
      {
         return m_v2;
      }
      //: Return feature row

      inline RealT Col2() const
      {
         return m_u2;
      }
      //: Return feature column
*/
      DListC< Tuple2C<IntT, Vector3dC> > AFeaturePoints(
         const SArray1dC<ImageC<ByteT> >& silimages,
         const MeshVisibilityC& vis,
         const IntT maxdist = 3,
         const RealT angthresh = 0.1) const;
      // Extract the feature points on the mesh corresponding to edge pixels in the input images
      // Returns the 3D position on the mesh and the index for the corresponding tri
      // maxdist defines the maximum distance from an edge pixel for a potential feature vertex
      // angthresh defines a threshold on the dot product between viewing direction and surface normal
      //  this helps to ensure feature points are caused by appearance rather than geometric discontinuities

      DListC< Tuple2C<IntT, Vector3dC> > AFeaturePoints2(
         const SArray1dC<ImageC<ByteT> >& silimages,
         const MeshVisibilityC& vis,
         const IntT camID,
         const IntT maxdist = 3,
         const RealT angthresh = 0.1) const;
      // Only for a single camera
      // Extract the feature points on the mesh corresponding to edge pixels in the input images
      // Returns the 3D position on the mesh and the index for the corresponding tri
      // maxdist defines the maximum distance from an edge pixel for a potential feature vertex
      // angthresh defines a threshold on the dot product between viewing direction and surface normal
      //  this helps to ensure feature points are caused by appearance rather than geometric discontinuities
      
      DListC< Tuple2C<IntT, Vector3dC> > AFeaturePoints3(
         const DListC<AFeature2C>& lp,
         const MeshVisibilityC& vis,
         const IntT camID) const;
      // Only for a single camera
      // Extract the feature points on the mesh corresponding to edge pixels in the input images
      // Returns the 3D position on the mesh and the index for the corresponding tri
      // maxdist defines the maximum distance from an edge pixel for a potential feature vertex
      // angthresh defines a threshold on the dot product between viewing direction and surface normal
      //  this helps to ensure feature points are caused by appearance rather than geometric discontinuities
/*
      SArray1dC<IntT> FeatVerts(
         const IntT camID,
         const DListC<AFeature2C>& lp,
         const SArray1dC<ImageC<ByteT> >& silimages,
         const MeshVisibilityC& vis,
         const IntT threshold = 50,
         const RealT angthresh = 0.7,
         const IntT edgedist = 3,
         const IntT maxdist = 3,
         const IntT nsize = 2,
         const IntT filter = 0) const;
      // Extract the feature vertices for the mesh
      // Corners are first extracted in the images
      // The corresponding visible mesh vertices are then found
      // edgedist defines the minimum distance from the silhouette edge for a valid edge pixel in the images
      // maxdist defines the maximum distance from an edge pixel for a potential feature vertex
      // angthresh defines a threshold on the dot product between viewing direction and surface normal
      //  this helps to ensure feature points are caused by appearance rather than geometric discontinuities
      // nsize sets the region size for a unique feature vertex, this extracts a single feature in a nsize neighbourhood
      //  this helps to ensure we do not get multiple corners at adjacent verts giving an ambiguous correspondence point
      // filter defines the number of times to smooth the image to extract corners at material rather than shading discontinuities

      //inline DListC< Tuple2C<IntT, Vector3dC> >& KPoints(
      //   const SArray1dC<ImageC<ByteT> >& silimages,
      //   const MeshVisibilityC& vis,
      //   const IntT& camID)
      //{
      //   return AFeaturePoints2(silimages,vis,camID);
      //}
      //: Return 3d feature points for single camera view
*/      
      void Load(istream &is, Index2dC topLeft = Index2dC(0,0));
      //: Load from stream.
   
//      void Load2(istream &is, Index2dC topLeft = Index2dC(0,0));
      //: Load from stream.

      void Save(ostream &os, Index2dC topLeft = Index2dC(0,0)) const;
      //: Save to stream.

   protected:

      static DListC< Tuple2C<IntT, Vector3dC> > AExtractMeshPoints(
         IntT camid,
         const ImageC<ByteT>& img,
         const PinholeCameraC& cam,
         const TriMeshC& mesh,
         const MeshVisibilityC& vis,
         const IntT maxdist,
         const RealT angthresh);
      // Determine the position on the mesh corresponding to a set of feature pixels in an image
      // Feature pixels are marked as 255 and all other pixels as 0
/*
      static DListC< Tuple2C<IntT, Vector3dC> > AExtractMeshPoints2(
         IntT camid,
         const DListC<AFeature2C>& lp,
         const PinholeCameraC& cam,
         const TriMeshC& mesh,
         const MeshVisibilityC& vis);
      // Returns 3d feature points from 2d locations for individual camera views
*/
      static DListC< Tuple2C<IntT, Vector3dC> > AExtractMeshPoints3(
         IntT camid,
	 const DListC<AFeature2C>& lp,
	 const PinholeCameraC& cam,
	 const TriMeshC& mesh,
	 const MeshVisibilityC& vis);
      // Returns 3d feature points from 2d locations for individual camera views

   protected:
      RealT m_u;
      RealT m_v;
//      RealT m_u1;
//      RealT m_v1;
//      RealT m_u2;
//      RealT m_v2;
      TriMeshC m_mesh;
      PinholeCameraArrayC m_cameras;
   };

   /*ostream &operator<<(ostream & s, const AFeature2C & p) {
      p.Save(s);
      return s;
   }*/
}

#endif 
