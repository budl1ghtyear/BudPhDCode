// Modified on 26/01/2009 by Ashish Doshi.
// Ravl includes

#include "Ravl/IO.hh"
#include "Ravl/DP/FileFormatIO.hh"

#include "Ravl/HSet.hh"
#include "Ravl/Polygon2d.hh"
#include "Ravl/Array2dIter.hh"
#include "Ravl/Array2dIter2.hh"
#include "Ravl/SArray1dIter2.hh"
#include "Ravl/Image/EdgeDetector.hh"
#include "Ravl/Image/CornerDetectorSusan.hh"
#include "Ravl/Image/DrawFrame.hh"
#include "Ravl/Image/BilinearInterpolation.hh"
#include "Ravl/3D/PinholeCamera0.hh"
#include "Ravl/3D/PinholeCamera1.hh"
#include "Ravl/3D/PinholeCamera2.hh"
#include "Ravl/3D/PinholeCameraArray.hh"

#include "Ravl/Image/Image.hh"
#include "Ravl/Image/ByteRGBValue.hh"
#include "Ravl/Image/PNMFormat.hh"

#include "../MultiViewToolkit/Matte/DistanceMap.hh"
using namespace MVTMatteN;

#include "../MultiViewToolkit/Texture/Texture.hh"
#include "../MultiViewToolkit/Texture/DepthImage.hh"
using namespace MVTTextureN;

#include "../MeshToolkit/MeshProc/MeshProc.hh"
using namespace MTMeshProcN;
#include "./AFeature2.hh"


// Use external prototypes
extern void InitPNGFormat();
extern void InitJPEGFormat();



namespace vbapAT
{

   AFeature2C::AFeature2C()
      : m_u(0.0),
        m_v(0.0)
   {
   }
   //: Default constructor

   AFeature2C::AFeature2C(RealT u,
                          RealT v)
      : m_u(u),
        m_v(v)
   {
   }
   //: Constructor

//   AFeature2C::AFeature2C(RealT u1,
//                          RealT v1,
//			  RealT u2,
//			  RealT v2)
//      : m_u1(u1),
//        m_v1(v1),
//        m_u2(u2),
//        m_v2(v2)
//   {
//   }
   //: Constructor

   // Determine the position on the mesh corresponding to a set of feature pixels in an image
   // Feature pixels are marked as 255 and all other pixels as 0
   DListC< Tuple2C<IntT, Vector3dC> >
      AFeature2C::AExtractMeshPoints(
      IntT camID,
      const ImageC<ByteT>& img,
      const PinholeCameraC& cam,
      const TriMeshC& mesh,
      const MeshVisibilityC& vis,
      const IntT maxdist,
      const RealT angthresh)
   {
      DListC< Tuple2C<IntT, Vector3dC> > ret;
      

      // Get the origin of the camera
      Vector3dC org;
      cam.Origin(org);

      // Create a distance map to the feature pixels
      ImageC<FloatT> dist = DistanceMap<ByteT>(img, (ByteT)255);
      //RavlN::Save("./testnew.png", dist);

      // Check each visible face to see if it contains a feature pixel
      IntT f = 0;
      IntT numfaces = mesh.Faces().Size();
      for (f = 0; f < numfaces; f++)
      {
         if (vis.FacetVisible(f, camID) || vis.FacetVisible(f, camID+1))
//         if (vis.FacetVisible(f, camID))
         {
            const TriC& tri = mesh.Faces()[f];
//            cerr << "Yes - 1" << endl;

            // Project the vertices onto the image
            // Check that a vertex on tri is within the distance threshold to a feature pixel
            bool bClose = false;
            DListC<Point2dC> points;
            for (IntT i = 0; i<3; i++)
            {
               Vector2dC p;
               cam.Project(p, tri.Position(i));
               points += Vector2dC(p.Y(), p.X());
               Index2dC pix(Round(p.Y()), Round(p.X()));
               if (dist.Frame().Contains(pix) && (IntT)dist[pix] <= maxdist)
               {
                  bClose = true;
//                  cerr << "Yes - 2" << endl;
               }

            }
            if (!bClose) continue;
            

            // Rasterise the tri and check for a feature pixel
            Polygon2dC poly(points);
            IndexRange2dC bbox = poly.BoundingRectangle().IndexRange();
//            cerr << img.Rows() << " " << img.Cols() << endl;
//            cerr << bbox.Rows() << " " << bbox.Cols() << endl;
            bbox.ClipBy(img.Frame()); // Clip by image size.
            // For each pixel inside bounding box...
            Array2dIterC<ByteT> it(img,bbox);
            IndexRange2dIterC itpix(bbox);
            for (; it; it++, itpix++) 
            {
//               cerr << "Yes - 3" << endl;

               // Check for a feature pixel
               if (it.Data())
               {
//                  cerr << "Yes - 4" << endl;
                  // Check if pixel is inside the polygon
                  Point2dC pnt(itpix.Data());
                  if (poly.Contains(pnt)) 
                  {
                     // Perform a plane - ray intersection to determine the position on the mesh
                     Vector3dC d;
                     cam.ProjectInverseDirection(d, Vector2dC(pnt.Y(), pnt.X()));
                     d.MakeUnit();
                     Vector3dC n = tri.FaceNormal();
                     //n *= -1; // Only use this if the triangles are the wrong way round.
                     if ((d.Dot(n)*-1.0) > angthresh)
                     {
                        RealT lambda = (n.Dot(tri.Position(0) - org)) / (n.Dot(d));
                        Vector3dC p = org + d * lambda;
                        ret += Tuple2C<IntT, Vector3dC>(f, p);
//                        cerr << mesh.Index(mesh.Faces()[f],0) << " " << mesh.Index(mesh.Faces()[f],1) << " " << mesh.Index(mesh.Faces()[f],2) << endl;
//                        cerr << tri[0] << tri[1] << tri[2] << endl;
//                        cerr << "Yes - 5" << endl;
//                        cerr << pnt.Y() << " " << pnt.X() << endl;
                     }
		     
                  }
               }
            }

         } // Visibility check
      } // foreach facet
      
      	
      return ret;
   }

/*
   // Determine the position on the mesh corresponding to a set of feature pixels for individual camera
   DListC< Tuple2C<IntT, Vector3dC> >
      AFeature2C::AExtractMeshPoints2(
      IntT camID,
      const DListC<AFeature2C>& lp,
      const PinholeCameraC& cam,
      const TriMeshC& mesh,
      const MeshVisibilityC& vis)
   {
             
      DListC< Tuple2C<IntT, Vector3dC> > ret;      
      SampleBilinearC<RealT, RealT> sampler;
      //const SArray1dC<TriC>& tris = mesh.Faces();
      //IntT numfaces = tris.Size();      

      // Get the origin of the camera
      Vector3dC org;
      cam.Origin(org);

      // Get depth information
      DepthImageC zrenderer(cam);
      zrenderer.Render(mesh);
      const ImageC<RealT>& zbuffer = zrenderer.Buffer();     

      {
         // Check for each feature point
         Vector3dC d;
         for(DLIterC<AFeature2C> lfit(lp); lfit; lfit++)
         {
            // Perform a plane - ray intersection to determine the position on the mesh
            Vector2dC imgPt(lfit->Col(), lfit->Row());
            Index2dC imgPix(Round(imgPt.Y()),Round(imgPt.X()));

            RealT depth;
            sampler(zbuffer, Vector2dC(imgPt.Y(), imgPt.X()), depth);
            cam.ProjectInverseDirection(d, Vector2dC(imgPt.Y(), imgPt.X()));
            d.MakeUnit();
            Vector3dC p = org + d*depth;
            //if (p.SqrEuclidDistance(s.Position()) < 0.0025)
            {
               ret += Tuple2C<IntT, Vector3dC>(lfit,p);
            }
         
         }
         

       }  
      	
      return ret;
   }
   // This particular function is NOT WORKING at the moment
*/
   
   // Determine the position on the mesh corresponding to a set of feature pixels for individual camera
   DListC< Tuple2C<IntT, Vector3dC> >
      AFeature2C::AExtractMeshPoints3(
      IntT camID,
      const DListC<AFeature2C>& lp,
      const PinholeCameraC& cam,
      const TriMeshC& mesh,
      const MeshVisibilityC& vis)
   {
      DListC< Tuple2C<IntT, Vector3dC> > ret;      
      SampleBilinearC<RealT, RealT> sampler;

      // Get the origin of the camera
      Vector3dC org;
      cam.Origin(org);

 
      // Check for each feature point
      Vector3dC d;
      for(DLIterC<AFeature2C> lfit(lp); lfit; lfit++)
      {
            Vector2dC imgPt(lfit->Col(), lfit->Row());
            Index2dC imgPix(Round(imgPt.Y()),Round(imgPt.X()));

            // Perform a plane - ray intersection to determine the position on the mesh
            Vector3dC d;
            RealT depth = 1;
            cam.ProjectInverseDirection(d, Vector2dC(imgPt.Y(), imgPt.X()));
            d.MakeUnit();
            {
               Vector3dC p = org + d*depth;
               ret += Tuple2C<IntT, Vector3dC>(lfit, p);
            }
		     
       }
         
       return ret;
   }
   // Converts the points from 2D to 3D - but not correctly to the point on the mesh

/*
   // Determine the position on the mesh corresponding to a set of feature pixels for individual camera
   DListC< Tuple2C<IntT, Vector3dC> >
      AFeature2C::AExtractMeshPoints4(
      IntT camID,
      const DListC<AFeature2C>& lp,
      const PinholeCameraC& cam,
      const TriMeshC& mesh,
      const MeshVisibilityC& vis)
   {
      DListC< Tuple2C<IntT, Vector3dC> > ret;      
      SampleBilinearC<RealT, RealT> sampler;

      // Get the origin of the camera
      Vector3dC org;
      cam.Origin(org);

 
      // Check for each feature point
      Vector3dC d;
      for(DLIterC<AFeature2C> lfit(lp); lfit; lfit++)
      {
            Vector2dC imgPt(lfit->Col(), lfit->Row());
            Index2dC imgPix(Round(imgPt.Y()),Round(imgPt.X()));

            // Perform a plane - ray intersection to determine the position on the mesh
            Vector3dC d;
            RealT depth = 4;
            cam.ProjectInverseDirection(d, Vector2dC(imgPt.Y(), imgPt.X()));
            d.MakeUnit();
            {
               Vector3dC p = org + d*depth;
               ret += Tuple2C<IntT, Vector3dC>(lfit, p);
            }
		     
       }
         
       return ret;
   }
   // Converts the points from 2D to 3D
*/

/*
   // Extract the feature vertices for the mesh
   // Returns the vertices closest to extracted feature pixels in the input images
   SArray1dC<IntT> 
      AFeature2C::FeatVerts(
         const IntT camID,
         const DListC<AFeature2C>& lp,
         const SArray1dC<ImageC<ByteT> >& silimages,
         const MeshVisibilityC& vis,
         const IntT threshold,
         const RealT angthresh,
         const IntT edgedist,
         const IntT maxdist,
         const IntT nsize,
         const IntT filter) const
   {
      HashC<IntT, RealT> fverts;
      const SArray1dC<TriC>& tris = m_mesh.Faces();
      IntT vsize = m_mesh.Vertices().Size();
      const ImageC<ByteT>& sdimage = silimages[camID]; 
//      DListC< Tuple2C<IntT, Vector3dC> > mpoints;

      // Find the vertices that are closest to edge pixels in the images
      //IntT cam;
      //IntT numimages = silimages.Size();
      //for (cam = 0; cam < numimages; cam++)
      {
         Vector3dC o;
         m_cameras[camID].Origin(o);
         //ImageC<ByteT> eimg = ExtractCornerImage(images[cam], silimages[cam], edgedist, filter, threshold);
         DListC< Tuple2C<IntT, Vector3dC> > fpts = AExtractMeshPoints(camID, silimages[camID], m_cameras[camID], m_mesh, vis, maxdist, angthresh);
//         DListC< Tuple2C<IntT, Vector3dC> > fpts = AExtractMeshPoints3(camID, lp, m_cameras[camID], m_mesh, vis);

         SArray1dC<IntT> fvts;
         SArray1dC<Vector3dC> mps = SArray1dC< Vector3dC >(lp.Size());
         SArray1dIterC< Vector3dC > itm(mps);
         SArray1dIterC< IntT > itv(fvts);         
         for(DLIterC<AFeature2C> lfit(lp); lfit; lfit++, itv++, itm++)
         {
            Vector2dC imgPt(lfit->Col(), lfit->Row());
//            Index2dC imgPix(Round(imgPt.Y()),Round(imgPt.X()));
        
            RealT eval = 200.0;
            IntT fv;
            Vector3dC pp;           
            DLIterC< Tuple2C<IntT, Vector3dC> > itp(fpts);            
            for (; itp; itp++)
            {
               Vector2dC featPt;
               m_cameras[camID].Project(featPt, itp.Data().Data2());
               Index2dC spix((IntT)(featPt.Y()), (IntT)(featPt.X()));
               if (!sdimage.Frame().Contains(spix))
               {
                  continue;
               }
               if (!(imgPt.EuclidDistance(featPt)<eval))
               {
                  continue;
               }
               else
               {
                  eval = imgPt.EuclidDistance(featPt);
                  pp = itp.Data().Data2();
//                cerr << itm.Data().X() << " " << itm.Data().Y() << " " << itm.Data().Z() << "\n"; 
                  fv = itp.Data().Data1();                     
               }               
            }
            itv.Data() = fv;
            itm.Data() = pp;
//            mpoints += Tuple2C<IntT, Vector3dC>(fv,pp);
         }
         cerr << mps.Size() << endl;
//         SArray1dIterC< Vector3dC > itmtest(mps);
//         for(; itmtest; itmtest++)
//         {
//            cerr << itmtest.Data().X() << " " << itmtest.Data().Y() << " " << itmtest.Data().Z() << "\n";
//         }
         
//         DListC<Tuple2C<IntT, Vector3dC> > mpoints(fvts, mps);
         
         // Iterate over the mesh points and extract the closest vertex to each point
         DLIterC< Tuple2C<IntT, Vector3dC> > it(fpts);
//         DLIterC< Tuple2C<IntT, Vector3dC> > it(mpoints);
         for (; it; it++)
         {
            // Find the closest vertex to the feature point
            const TriC& tri = tris[it.Data().Data1()];
            RealT dist = it.Data().Data2().SqrEuclidDistance(tri.Position(0));
            IntT closest = 0;
            for (IntT i=1; i<3; i++)
            {
               RealT d = it.Data().Data2().SqrEuclidDistance(tri.Position(i));
               if (d < dist)
               {
                  dist = d;
                  closest = i;
               }
            }
            IntT v = m_mesh.Index(tri, closest); ;
//            valid += v;
            Vector3dC dir = Vector3dC(o - it.Data().Data2()).Unit();
            RealT nval = m_mesh.Vertices()[v].Normal().Dot(dir);
            RealT norm;
            if (fverts.Lookup(v, norm))
            {
               if (norm < nval)
               {
                  fverts[v] = nval;
               }
            }
            else
            {
               fverts[v] = nval;
            }
         }
      }

      // Now decimate vertices with a closer normal in the nsize neighbourhood
      DListC<IntT> valid;
      SArray1dC< DListC<UIntT> > vneighbourhood = VertexNeighbourhood(m_mesh);
      SArray1dC<IntT> vdist(vsize);
      HashIterC<IntT,RealT> itfv(fverts);
      for (; itfv; itfv++)
      {
         IntT v = itfv.Key();
         RealT norm = itfv.Data();
         // Check the nsize neighbourhood
         bool bValid = true;
         vdist.Fill(nsize+1);
         DListC< IntT > queue;
         queue += v;
         vdist[v] = 0;
         DLIterC< IntT > itq(queue);
         for (; itq && bValid; itq++)
         {
            IntT vn = itq.Data();
            IntT ndepth = vdist[vn] + 1;
            // Visit the 1-neighbourhood of this vertex
            DLIterC<UIntT> itn( vneighbourhood[vn] );
            for (; itn; itn++)
            {
               IntT nid = itn.Data();
               // Check whether the neighbour has been visited already
               if (vdist[nid] > ndepth)
               {
                  RealT nnorm;
                  if (fverts.Lookup(nid, nnorm))
                  {
                     if (nnorm > norm)
                     {
                        bValid = false;
                        break;
                     }
                  }
                  queue += nid;
                  vdist[nid] = ndepth;
               }
            }
         }
         if (bValid)
         {
            valid += v;
         }
      }

      IntT vcount = valid.Size();
      SArray1dC<IntT> ret(vcount);
      SArray1dIterC<IntT> it(ret);
      DLIterC<IntT> itv(valid);
      for (; itv; itv++, it++)
      {
         it.Data() = itv.Data();
      }

      return ret;
   }
*/


   // Extract the feature points on the mesh corresponding to pixels in the input images
   // Returns the 3D position on the mesh and the index for the corresponding tri
   DListC< Tuple2C<IntT, Vector3dC> >
      AFeature2C::AFeaturePoints(
         const SArray1dC<ImageC<ByteT> >& silimages,
         const MeshVisibilityC& vis,
         const IntT maxdist,
         const RealT angthresh) const
   {
      DListC< Tuple2C<IntT, Vector3dC> > ret;

      // Extract the mesh points for each edge pixel across all views
      IntT cam;
      IntT numimages = silimages.Size();
      for (cam = 0; cam < numimages; cam++)
      {
         ret += AExtractMeshPoints(cam, silimages[cam], m_cameras[cam], m_mesh, vis, maxdist, angthresh);
      }
      
      return ret;
   }
   //: Load 3d feature points

   // Extract the feature points on the mesh corresponding to pixels in the input images
   // Returns the 3D position on the mesh and the index for the corresponding tri
   DListC< Tuple2C<IntT, Vector3dC> >
      AFeature2C::AFeaturePoints2(
         const SArray1dC<ImageC<ByteT> >& silimages,
         const MeshVisibilityC& vis,
         const IntT camID,
         const IntT maxdist,
         const RealT angthresh) const
   {
      DListC< Tuple2C<IntT, Vector3dC> > ret;

      // Extract the mesh points for each edge pixel across all views
      ret = AExtractMeshPoints(camID, silimages[camID], m_cameras[camID], m_mesh, vis, maxdist, angthresh);      
      return ret;
   }
   //: Load 3d feature points for single camera

   // Returns the 3D positions on the mesh from 2D feature locations in an image.
   DListC< Tuple2C<IntT, Vector3dC> >
      AFeature2C::AFeaturePoints3(        
	 const DListC<AFeature2C>& lp,
         const MeshVisibilityC& vis,
         const IntT camID) const
   {
      DListC< Tuple2C<IntT, Vector3dC> > ret;

      // Extract the mesh points for each edge pixel across all views
      ret = AExtractMeshPoints3(camID, lp, m_cameras[camID], m_mesh, vis);      
      return ret;
   }
   //: Load 3d feature points for single camera
/*
   // Returns the 3D positions on the mesh from 2D feature locations in an image.
   DListC< Tuple2C<IntT, Vector3dC> >
      AFeature2C::AFeaturePoints4(        
	 const DListC<AFeature2C>& lp,
         const MeshVisibilityC& vis,
         const IntT camID) const
   {
      DListC< Tuple2C<IntT, Vector3dC> > ret;
      DListC< Tuple2C<IntT, Vector3dC> > ret2;
      const SArray1dC<TriC>& tris = m_mesh.Faces();
      IntT vsize = m_mesh.Vertices().Size();

      // Extract the mesh points for each edge pixel across all views
      ret = AExtractMeshPoints3(camID, lp, m_cameras[camID], m_mesh, vis);

      DLIterC< Tuple2C<IntT, Vector3dC> it(ret);
      for (; it; it++)
      {
         const TriC& tri = tris[it.Data().Data1()];
         RealT dist = it.Data().Data2().SqrEuclidDistance(tri.Position(0));
         IntT closest = 0;
         
       
      return ret;
   }
   //: Load 3d feature points for single camera
*/
   void AFeature2C::Load(istream &is, Index2dC topLeft)
   {
      is >> m_u >> m_v;
      m_u += topLeft.Col();
      m_v += topLeft.Row();
   }
   //: Load feature locations

//   void AFeature2C::Load2(istream &is, Index2dC topLeft)
//   {
//      is >> m_u1 >> m_v1 >> m_u2 >> m_v2;
//      m_u1 += topLeft.Col1();
//      m_v1 += topLeft.Row1();
//      m_u2 += topLeft.Col2();
//      m_v2 += topLeft.Row2();
//   }
   //: Load feature locations

   void AFeature2C::Save(ostream &os, Index2dC topLeft) const
   {
      os << m_u-topLeft.Col() <<  " " << m_v-topLeft.Row();
   }
}
