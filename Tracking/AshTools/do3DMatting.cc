// RAVL includes
#include "Ravl/IO.hh"
#include "Ravl/DP/FileFormatIO.hh"
#include "Ravl/Hash.hh"
#include "Ravl/Option.hh"
#include "Ravl/Tuple2.hh"
#include "Ravl/Point2d.hh"
#include "Ravl/Array1d.hh"
#include "Ravl/HashIter.hh"
#include "Ravl/Array1dIter.hh"
#include "Ravl/Array1dIter2.hh"
#include "Ravl/Array1dIter3.hh"
#include "Ravl/Array1dIter4.hh"
#include "Ravl/Array1dIter5.hh"
#include "Ravl/Array2dIter2.hh"
#include "Ravl/Array2dIter3.hh"
#include "Ravl/Array2dIter4.hh"
#include "Ravl/Array2dIter5.hh"
#include "Ravl/Image/ImgIO.hh"
#include "Ravl/3D/TriMesh.hh"
#include "Ravl/3D/RavlMeshIO.hh"
#include "Ravl/3D/PinholeCameraArray.hh"
#include "Ravl/Index3d.hh"
#include "Ravl/DList.hh"
#include "Ravl/DLIter.hh"
#include "Ravl/Image/WarpScale.hh"
#include "Ravl/Image/DrawCross.hh"
#include "Ravl/Image/DrawCircle.hh"
#include "Ravl/TFVector.hh"
#include "Ravl/LeastSquares.hh"
#include "Ravl/Image/GaussConvolve2d.hh"
#include "Ravl/Image/RealRGBValue.hh"
#include "Ravl/Image/ByteRGBValue.hh"
#include "Ravl/Image/ByteRGBAValue.hh"
#include "Ravl/Image/SqrComposition.hh"
#include "Ravl/Image/EdgeDeriche.hh"
#include "Ravl/Image/EdgeSobel.hh"
#include "Ravl/Image/EdgeNonMaxSuppression.hh"
#include "Ravl/Image/Edgel.hh"
#include "Ravl/Image/EdgeLink.hh"
#include "Ravl/Image/Dilate.hh"
#include "Ravl/Image/Erode.hh"
#include "Ravl/Assert.hh"
#include "Ravl/Sums1d2.hh"
#include "Ravl/SumsNd2.hh"
#include "Ravl/Collection.hh"

using namespace RavlN;
using namespace Ravl3DN;
using namespace RavlImageN;
using namespace RavlConstN;

// Local includes
#include "MVT/Matte/Matte.hh"
#include "MVT/Matte/BayesianMatte.hh"
using namespace MVTMatteN;
#include "MVT/Stereo/Stereo.hh"
//#include "MVT/Stereo/StereoPair.hh"
using namespace MVTStereoN;
#include "MVT/Texture/DepthImage.hh"
using namespace MVTTextureN;

// MTMeshProc includes
#include "MT/MeshProc/MeshProc.hh"
#include "MT/MeshProc/MeshSmooth.hh"
using namespace MTMeshProcN;

#include "MaxFlow/graph.h"

#include "MultiViewMatching.hh"
#include "BinarySkeletonisation.hh"
#include "BinaryPruning.hh"

#include "Compute2DMatte.hh"


typedef TFVectorC< PairC<IntT>,5> RangeVectorC;
const IntT invalidIndex = 0;
const ByteRGBValueC BGPixel(0,0,0);
const RealT maxScore = 100.0;

class OrderedVectorC
  : public VectorC
{
  public:
    OrderedVectorC()
      : VectorC()
    {}
    //: Default constructor.
    
    explicit OrderedVectorC(SizeT size)
      : VectorC(size)
    {}
    //: Constructor.
    // Create a vector of 'size' elements
    
    OrderedVectorC(const VectorC &oth)
      : VectorC(oth)
    {}
    //: Base class constructor.

    bool operator<(const OrderedVectorC &v) 
    { return(SumOfAbs() < v.SumOfAbs()); }
    //: Comparison operator

    bool operator>(const OrderedVectorC &v) 
    { return(SumOfAbs() > v.SumOfAbs()); }
    //: Comparison operator
};

bool RobustMeanCovarianceC(const DListC<VectorC> & data,
                           VectorC& mean,
                           MatrixRSC& cov,
                           bool zeroMean);

bool
ComputeEdgeMap(const ImageC<ByteRGBValueC>& rgbImage,
               ImageC<bool>& edgeMap,
               RealT minHyst = 15,
               RealT maxHyst = 20,
               bool useSobel = false,
               RealT omegaDeriche = 0.001,
               RealT alphaDeriche = 2.0);
// Compute edge map

bool
ComputeProjections(SArray1dC<Vector2dC>& imgPts,
                   const Vector3dC& Pt,
                   const PinholeCameraArrayC& cameras,
                   const SArray1dC<IntT>& camIDs);
// Determine the projection of the point Pt into the images with IDs CamIDs

bool
ComputeProjections(SArray1dC<Vector2dC>& imgPts,
                   DListC<IntT>& fgCamIDs,
                   DListC<IntT>& bgCamIDs,
                   DListC<IntT>& ocCamIDs,
                   const Vector3dC& Pt,
                   const PinholeCameraArrayC& cameras,
                   const SArray1dC<IntT>& camIDs,
                   const SArray1dC<ImageC<ByteT> >& triMaps,
                   const IndexRange2dC activeWindow);
// Determine the projection of the point Pt into the images with IDs CamIDs,
//  partitioning cameras into foreground and background

bool
IsWithinVisualHull(const Vector3dC& Pt,
                   const PinholeCameraArrayC& cameras,
                   const SArray1dC<ImageC<ByteT> >& triMaps,
                   const IndexRange2dC activeWindow,
                   bool verbose = false);
// Determines if a 3D point is within the visual hull by looking at
//  its projection in the visual hull silhouettes

bool
ComputeIntersectionRange(ImageC<RangeVectorC>& VHrange,
                         ImageC<PairC<IntT> >& BGrange,
                         const SArray1dC<ImageC<ByteT> >& triMaps,
                         const ImageC<bool>& refTriMapMask,
                         const SArray1dC<ImageC<RealT> >& depthImgs,
                         const PinholeCameraArrayC& cameras,
                         const IndexRange2dC activeWindow,
                         IntT refCamID,
                         RealT yPitchMin,
                         RealT yPitchMax,
                         RealT step,
                         IntT maxSurfDist = 1e9);
// Determine the intersection range for camera image

bool
ComputeVisibility(SArray1dC<IntT>& visCamIDs,
                  const SArray1dC<ImageC<RealT> >& depthImgs,
                  const Vector3dC& Pt,
                  const PinholeCameraArrayC& cameras,
                  const SArray1dC<IntT>& camIDs,
                  RealT tolVis);
// Determine the visibility of the point Pt in the images with IDs CamIDs

bool SmoothnessCost(RealT& score,
                    bool Edge1,
                    bool Edge2,
                    RealT kSmoothMin = 0.1,
                    RealT kSmoothMax = 1.0);
// Return smoothness cost adapted according to edge detection

bool ComputeFGScoreMean(RealT& score,
                   const ImageC<ByteRGBValueC>& refPatch,
                   const SArray1dC<ImageC<ByteRGBValueC> >& oPatch,
                   const ImageC<ByteRGBValueC>& refBGPatch,
                   const SArray1dC<ImageC<ByteRGBValueC> >& oBGPatch,
                   RealT tData,
                   RealT kData,
                   SArray1dC<Vector3dC> fgVars,
                   Vector3dC refbgVar,
                   SArray1dC<Vector3dC> obgVar,
                   IntT maxBackground=(IntT)1e6,
                   bool VERBOSE = false);
// Return the sum of squarred difference between a reference image window and other image windows
//  using background information
// Returns false if the number of background pixels exceeds the tolerance

bool ComputeFGScoreBest(RealT& score,
                   const ImageC<ByteRGBValueC>& refPatch,
                   const SArray1dC<ImageC<ByteRGBValueC> >& oPatch,
                   const ImageC<ByteRGBValueC>& refBGPatch,
                   const SArray1dC<ImageC<ByteRGBValueC> >& oBGPatch,
                   RealT tData,
                   RealT kData,
                   SArray1dC<Vector3dC> fgVars,
                   Vector3dC refbgVar,
                   SArray1dC<Vector3dC> obgVars,
                   IntT maxBackground=(IntT)1e6,
                   bool VERBOSE = false);
// Return the sum of squarred difference between a reference image window and other image windows
//  using background information
// Returns false if the number of background pixels exceeds the tolerance

bool SSDCorrelation(RealT& score,
                    const ImageC<ByteRGBValueC>& keyPatch,
                    const ImageC<ByteRGBValueC>& offsetPatch,
                    const ImageC<ByteRGBValueC>& keyBGPatch,
                    const ImageC<ByteRGBValueC>& offsetBGPatch,
                    RealT tData,
                    RealT kData,
                    Vector3dC fgVar = Vector3dC(1.0,1.0,1.0),
                    PairC<Vector3dC> bgVars = PairC<Vector3dC>(Vector3dC(1.0,1.0,1.0), Vector3dC(1.0,1.0,1.0)),
                    IntT maxBackground=(IntT)1e6,
                    bool VERBOSE = false);
// Return the sum of squarred difference between two image windows
//  using background information
// Returns false if the number of background pixels exceeds the tolerance

bool SSDCorrelation(RealT& score,
                    const ImageC<ByteRGBValueC>& keyPatch,
                    const ImageC<ByteRGBValueC>& offsetPatch,
                    Vector3dC var = Vector3dC(1.0,1.0,1.0),
                    IntT maxBackground=(IntT)1e6,
                    bool VERBOSE = false);
// Return the sum of squarred difference between two image windows
// Returns false if the number of background pixels exceeds the tolerance

bool comp_le(const RealT & el1, const RealT& el2) {
  return(el1 <=  el2);
}


RealT
GraphCut(ImageC<RealT>& depthMap,
         const ImageC<RangeVectorC>& VHrange,
         const ImageC<PairC<IntT> >& BGrange,
         const SArray1dC<ImageC<ByteRGBValueC> >& rgbImages,
         const SArray1dC<ImageC<ByteRGBValueC> >& bgmeanImages,
         const SArray1dC<ImageC<RealT> >& depthImgs,
         const PinholeCameraArrayC& cameras,
         IntT refCamID,
         const SArray1dC<IntT>& oCamIDs,
         const SArray1dC<Vector3dC>& fgVars,
         const SArray1dC<ImageC<RealRGBValueC> >& bgVarImages,
         const ImageC<bool>& edgeMap,
         RealT step,
         RealT tolVis,
         bool useMean,
         IntT winsize,
         RealT kaSmoothMin,
         RealT kaSmoothMax,
         RealT kdSmoothMin,
         RealT kdSmoothMax,
         RealT atData,
         RealT sigmaFG = 1.0,
         RealT sigmaBG = 1.0,
         RealT kData = 5.0,
         IntT maxBackground = (IntT)1e6,
         bool noBG = false,
         IndexRange2dC activeWindow = IndexRange2dC(),
         bool useBestScore=false,
         bool FG_VERBOSE = false,
         bool BG_VERBOSE = false,
         Index2dC testImgPix = Index2dC(-1e6, -1e6));
// Perform a graph-cut to maximise photoconsistency between views
//  imposing smoothness constaint based on edginess

bool
BuildMesh(const ImageC<RealT>& depthMap,
          const SArray1dC<ImageC<ByteT> >& triMaps,
          const PinholeCameraArrayC& cameras,
          const IndexRange2dC activeWindow,
          IntT refCamID,
          RealT step,
          IntT smooth,
          const ImageC<ByteRGBValueC>& fgRefImage,
          bool embedColour,
          IntT frameNumber,
          StringC output,
          const ImageC<ByteT>& alphaMatte = ImageC<ByteT>());
// Builds a mesh of the reconstruction
//  the mesh is split at regions where the depth discontinuity is larger than depthLimit

bool
AreConnected(const PairC<Vector3dC>& pts,
             const SArray1dC<ImageC<ByteT> >& triMaps,
             const PinholeCameraArrayC& cameras,
             const IndexRange2dC activeWindow,
             RealT step);

bool
InterpolateColor(const ImageC<ByteRGBValueC>& rgbImage,
                 const Vector2dC& imgPt,
                 RealRGBValueC& colour);
// Interpolate colour from four pixel neighbours

bool
ProjectToReference(const SArray1dC<ImageC<ByteRGBValueC> >& rgbImages,
                   const ImageC<RealT>& depthMap,
                   const PinholeCameraArrayC& cameras,
                   IntT refCamID,
                   SArray1dC<IntT>& auxCamIDs,
                   SArray1dC<ImageC<RealRGBValueC> >& warpedImages,
                   SArray1dC<ImageC<RealT> >& errorMaps,
                   const IndexRange2dC activeWindow,
                   IntT frameNumber,
                   StringC outputSTATS,
                   bool useRobustStats = false,
                   bool useIsotropic = false,
                   bool VERBOSE = false);
// Project auxiliary cameras to the reference camera

int main(int argc, char* argv[])
{
   // Link image formats
   InitDPImageIO();
   InitTriFormat();
   InitTriMeshIO();

   // Command line options ///////////////////////////////////////////////////

   // Get the command line
   OptionC option(argc,argv);

   // INPUT DATA PARAMETERS
   StringC inputCAM        = option.String("cal","data/calibration/0500.cal","Camera calibration file");
   StringC inputRGB        = option.String("rgb","data/%02d/originals/%04d.tif","Input images");
   StringC inputTRI        = option.String("tri","data/%02d/trimaps/%04d.tif","Input trimaps");
   StringC inputMESH       = option.String("m","reconstructions/cons_vh/%d.tri","Visual hull");
   StringC bginputSTATS    = option.String("bgstatsin","reconstructions/graph_cut/stats/%02d/%02d.%04d.txt","Input stats for bg");
   StringC fginputSTATS    = option.String("fgstatsin","reconstructions/graph_cut/stats/%02d/%04d.txt","Input stats for fg");
   RealT fgcor             = option.Real("fgcor", 1.0, "fg correction for variance");

   // OUTPUT DATA PARAMETERS
   StringC outputDEPTH     = option.String("dres","reconstructions/graph_cut/depth_maps/%02d/%04d.png","Output depth map of foreground reconstruction");
   StringC outputMESH      = option.String("mres","reconstructions/graph_cut/meshes/%02d/%04d.tri","Output mesh of foreground reconstruction");
   StringC outputDEPTH2    = option.String("dres2","reconstructions/graph_cut/depth_maps_opt/%02d/%04d.png","Output depth map of foreground reconstruction (after boundary refinement)");
   StringC outputMESH2     = option.String("mres2","reconstructions/graph_cut/meshes_opt/%02d/%04d.tri","Output mesh of foreground reconstruction (after boundary refinement)");
   StringC outputMATTE     = option.String("ares","reconstructions/graph_cut/alpha_mattes/%02d/%04d.tri","Output alpha matte");
   StringC outputFG        = option.String("fgres","reconstructions/graph_cut/foreground/%02d/%04d.tri","Output foreground image");
   StringC outputSTATS     = option.String("statsres","","Output statistics");

   // MAIN PARAMETERS
   StringC inputIDs        = option.String("ids","0,6,7","List of IDs of the cameras considered, the first camera will be used as the reference camera (e.g. 5,4,6,3,7)");
   IntT frameNumber        = option.Int("f", 500, "frame number");
   RealT depthStep         = option.Real("dstep", 0.05, "Sampling step for depth");
   IntT smooth             = option.Boolean("smooth",0,"Mesh smoothing (number of iterations)");
   bool VH                 = option.Boolean("vh",false,"Use visual hull to compute depth");
   bool VH2                = option.Boolean("vh2",false,"Use visual hull to compute depth");
   Index2dC org            = option.Index2d("roiorg", 0, 0, "Origin of region of interest in image");
   Index2dC end            = option.Index2d("roiend", -1, -1, "End of region of interest in image");
   RealT sigmaFG           = option.Real("sigmafg", 2.0, "Standard deviation of matched foreground colours");
   RealT sigmaBG           = option.Real("sigmabg", 1.0, "Standard deviation of background colours");
   IntT maxSurfDist        = option.Int("msd", (IntT)1e6, "Maximum distance from visual hull surface (in depth steps)");
   bool REFINE_BOUNDARY    = option.Boolean("refine_boundary",false,"Perform matte refinement on player boundary strips");
   IntT cuttingVal         = option.Int("cut", 0, "Minimum alpha value allowed");
   bool embedColour        = option.Boolean("embedColour",false,"Save output as a coloured mesh");
   bool definiteFGonly     = option.Boolean("definitefgonly",false,"Reconstruct only definite foreground");

   // BACKGROUND PARAMETERS
   StringC inputBGmean     = option.String("bgm","data/%02d/background.tif","Input background mean images");
   StringC inputBGstd      = option.String("bgv","","Input background standard deviation images");
   bool noBG               = option.Boolean("nobg",false,"Do not include background layer");
   RealT yPitchMin         = option.Real("ymin", -0.1, "Minimum pitch altitude");
   RealT yPitchMax         = option.Real("ymax", 0.0, "Maximum pitch altitude");

   // DATA COST RELATED PARAMETERS
   IntT winsize            = option.Int("win", 1, "Window size e.g 2 = 5x5 window");
   IntT maxBackground      = option.Int("maxbg", Floor(0.8*Sqr(2*winsize+1)), "Maximum number of pixels in the window allowed to be background for matching");
   RealT tData             = option.Real("tdata", 0.001, "Threshold for matching with BG");
   RealT kData             = option.Real("kdata", 3.0, "BG penalty");
   RealT tolVis            = option.Real("tolvis", 0.25, "Tolerance on visibility constraint");
   bool noVis              = option.Boolean("novis",false,"Do not check visibility");
   bool useMean            = option.Boolean("mean",false,"Use mean for matching score");
   bool useBestScore       = option.Boolean("best",false,"Use best of separate foreground matching scores");

   // SMOOTHNESS COST RELATED PARAMETERS
   RealT kaSmoothMin       = option.Real("kasmin", 0.002, "minimum smoothing weight for alpha");
   RealT kaSmoothMax       = option.Real("kasmax", 0.004, "maximum smoothing weight for alpha");
   RealT kdSmoothMin       = option.Real("kdsmin", 0.002, "minimum smoothing weight for d");
   RealT kdSmoothMax       = option.Real("kdsmax", 0.003, "maximum smoothing weight for d");
   bool useSobel           = option.Boolean("sobel",false,"Use Sobel filter instead of Deriche filter");
   RealT omegaDeriche      = option.Real("oDeriche", 0.001, "omega parameter for Deriche filter");
   RealT alphaDeriche      = option.Real("aDeriche", 2.0, "alpha parameter for Deriche filter");
   RealT minHyst           = option.Real("minHyst", 5.0, "minimum hysteresis threshold for edge detection");
   RealT maxHyst           = option.Real("maxHyst", 10.0, "maximum hysteresis threshold for edge detection");
   IntT dilate             = option.Int("dilate",0,"Size of dilation kernel (n means a kernel of size 2*n+1) to apply to edges");

   IntT dilate_ref         = option.Int("dilate_ref",1,"Size of dilation kernel for refinement");
   IntT erode_ref          = option.Int("erode_ref",2,"Size of erosion kernel for refinement");
   IntT pruning_ref        = option.Int("pruning_ref",1,"Size of pruning kernel for refinement");

   // VERBOSITY PARAMETERS
   bool VERBOSE            = option.Boolean("v", false, "Verbose mode");
   bool VERBOSE_FG         = option.Boolean("vfg", false, "Verbose mode for foreground");
   bool VERBOSE_BG         = option.Boolean("vbg", false, "Verbose mode for background");
   Index2dC testImgPix     = option.Index2d("test", (IntT)-1e6, (IntT)-1e6, "Coordinate of image test point");

   bool useRobustStats     = option.Boolean("robust",false,"Compute robust estimate of mean and covariance matrix based on median");
   bool useIsotropic       = option.Boolean("iso",false,"Assume isotropic Gaussian distribution");

   option.Check();

   // Load Data //////////////////////////////////////////////////

   char pFile[255];

   // Load the camera calibration
   PinholeCameraArrayC cameras;
   IStreamC camFile(inputCAM);
   camFile >> cameras;
   camFile.Close();
   IntT numCams = cameras.Size();

   // Read list of camera IDs
   SArray1dC<IntT> allCamIDs = GetMultiViewIDs(inputIDs);

   // The first camera is used for reference (principal view)
   IntT refCamID = allCamIDs[0];
   SArray1dC<IntT> oCamIDs(allCamIDs, allCamIDs.Size()-1, 1);

   // Load stats
   SArray1dC<Vector3dC> fgVars(oCamIDs.Size());
   fgVars.Fill(Vector3dC(1.0,1.0,1.0));
   if(fginputSTATS.Size()>0) {
      // read foreground statistics for all pairs used in reconstuction
      for(Array1dIter2C<IntT, Vector3dC> it(oCamIDs, fgVars); it; it++) {
         Vector3dC fgvar;
         sprintf(pFile,fginputSTATS,refCamID,it.Data1(),frameNumber);
         IStreamC statsFile(pFile);
         statsFile >> fgvar;
         statsFile.Close();
         it.Data2() = fgcor*Sqr(fgvar/255.0);
      }
   }

   if(VERBOSE) {
      cerr << "fgVars: " << fgVars << endl;
   }

   // Load trimaps
   SArray1dC<ImageC<ByteT> > triMaps(numCams);

   // Load the visual hull mesh and generate trimaps automatically for all cams except ref cam
   TriMeshC mesh;
   sprintf(pFile,inputMESH,frameNumber);
   if(!Load(pFile, mesh, "", VERBOSE))
   {
      return 1;
   }
   if(noVis) {
      cerr << "WARNING: no visibility handling" << endl;
   }
   cerr << "Computing depth maps";
   SArray1dC<ImageC<RealT> > depthImgs(numCams);
   for (IntT cam=0; cam<numCams; cam++) 
   {
      // Render the depth map in each camera
      DepthImageC zbuffer(cameras[cam]);
      zbuffer.Render(mesh, false);
      depthImgs[cam] = zbuffer.Buffer();
      cerr << ".";

      triMaps[cam] = ImageC<ByteT>(depthImgs[cam].Frame());
      triMaps[cam].Fill(0);
      for(Array2dIter2C<RealT, ByteT> it(depthImgs[cam], triMaps[cam]); it; it++) {
         if(it.Data1()>0.0) {
            it.Data2() = 127;
         }
      }
   }
   cerr << endl;

   // load trimap for ref cam
   sprintf(pFile,inputTRI,refCamID,frameNumber);
   if(!Load(pFile,triMaps[refCamID],"",VERBOSE)) 
   {
      return 1;
   }

   // temporary added to test the effect of no having definite fg areas
   /*for (Array2dIterC<ByteT> it2(triMaps[refCamID]); it2; it2++)
   {
      if (it2.Data() ==255)
      {
         it2.Data() = 127;
      }
   }*/

   if(definiteFGonly) {

      noBG = true;

      for (Array2dIterC<ByteT> it2(triMaps[refCamID]); it2; it2++)
      {
         if (it2.Data() !=255)
         {
            it2.Data() = 0;
         }
      }
   }

   // Load the rest of the images needed for reconstruction
   SArray1dC<ImageC<ByteRGBValueC> > rgbImages(numCams); // input images
   SArray1dC<ImageC<ByteRGBValueC> > fgImages(numCams); // foreground images
   SArray1dC<ImageC<ByteRGBValueC> > bgmeanImages(numCams); // backgound images
   SArray1dC<ImageC<RealRGBValueC> > bgVarImages(numCams); // backgound images
   for (Array1dIterC<IntT> it(allCamIDs); it; it++) 
   {
      IntT cam = it.Data();

      sprintf(pFile,inputRGB.chars(),cam,frameNumber);
      if(!Load(pFile,rgbImages[cam],"",VERBOSE)) 
      {
         return 1;
      }
      if(!Load(pFile,fgImages[cam],"",VERBOSE)) 
      {
         return 1;
      }
      sprintf(pFile,inputBGmean.chars(),cam,frameNumber);
      if(!Load(pFile,bgmeanImages[cam],"",VERBOSE)) 
      {
         return 1;
      }
      if(!inputBGstd.IsEmpty()) {
         sprintf(pFile,inputBGstd.chars(),cam,frameNumber);
         if(!Load(pFile,bgVarImages[cam],"",VERBOSE)) 
         {
            return 1;
         }
      } else {
         Vector3dC bgstd;
         if(bginputSTATS.Size()>0) {
            sprintf(pFile,bginputSTATS.chars(),cam,frameNumber);
            IStreamC statsFile(pFile);
            statsFile >> bgstd;
            statsFile.Close();
         } else {
            bgstd = Vector3dC(1.0,1.0,1.0);
         }
         if(VERBOSE) {
            cerr << bgstd[0] << endl;
         }
         bgVarImages[cam] = ImageC<RealRGBValueC>(bgmeanImages[cam].Frame());
         bgVarImages[cam].Fill(RealRGBValueC(bgstd[0],bgstd[1],bgstd[2]));
      }
      // normalise and convert to variances
      for (Array2dIterC<RealRGBValueC> it2(bgVarImages[cam]); it2; it2++)
      {
         it2.Data() = Sqr(it2.Data()/255.0);
      }

      // Mask the RGB images to eliminate background so that we only perform stereo on the foreground
      for (Array2dIter2C<ByteRGBValueC, ByteT> it2(fgImages[cam], triMaps[cam]); it2; it2++)
      {
         if (it2.Data2() == 0)
         {
            it2.Data1() = BGPixel;
         }
      }
   }

   const ImageC<ByteRGBValueC>& rgbRefImage = rgbImages[refCamID];
   const ImageC<ByteRGBValueC>& fgRefImage = fgImages[refCamID];
   const ImageC<ByteRGBValueC>& bgmeanRefImage = bgmeanImages[refCamID];
   const ImageC<RealRGBValueC>& bgstdRefImage = bgVarImages[refCamID];
   const ImageC<ByteT>& refTriMap = triMaps[refCamID];

   // Dilate the non-background region to create a background layer around it and create a mask
   //  The background prior will be enforce on this layer instead of the full background
   //  This is mathematically equivalent for the graph-cut minimisation and more memory efficient with HD images
   ImageC<bool> refTriMapMask(refTriMap.Frame());
   ImageC<ByteT> bgkernel(IndexRange2dC(-2,2,-2,2));
   bgkernel.Fill(1);
   ImageC<ByteT> tmpTriMap(refTriMap.Frame());
   tmpTriMap.Fill(0);
   BinaryDilate(refTriMap.Copy(), bgkernel, tmpTriMap, ByteT(1));
   // Copy layer to trimap setting value to 1
   refTriMapMask.Fill(false);
   Array2dIter2C<ByteT, bool> it2(tmpTriMap, refTriMapMask);
   for (; it2; it2++)
   {
      if (it2.Data1()==1)
      {
         it2.Data2() = true;
      }
   }

   if(VERBOSE) {
      Save("@X:reference input image", rgbRefImage, "", VERBOSE);
      Save("@X:reference foreground image", fgRefImage, "", VERBOSE);
      Save("@X:reference mean background image", bgmeanRefImage, "", VERBOSE);
      Save("@X:reference std background image", bgstdRefImage, "", VERBOSE);
      Save("@X:reference trimap", refTriMap, "", VERBOSE);
   }

   ImageC<bool> edgeMap(rgbRefImage.Frame());
   ComputeEdgeMap(rgbRefImage, edgeMap, minHyst, maxHyst, useSobel, omegaDeriche, alphaDeriche);

   if(dilate>0) {
      // dilate edges
      ImageC<bool> kernel(IndexRange2dC(-dilate,dilate,-dilate,dilate));
      kernel.Fill(false);
      for(Array2dIterC<bool> it(kernel); it; it++) {
         if(Sqr(it.Index().Row())+Sqr(it.Index().Col()) <= Sqr(dilate)) {
            it.Data() = true;
         }
      }
      BinaryDilate(edgeMap.Copy(), kernel, edgeMap, true);
   }

   if(VERBOSE) {
      Save("@X:edge map", edgeMap, "", VERBOSE);

      // overlay edges in red on reference image
      ImageC<ByteRGBValueC> overlay(rgbRefImage.Copy());
      for(Array2dIter2C<ByteRGBValueC, bool> it(overlay, edgeMap); it; it++) {
         if(it.Data2()==true) {
            it.Data1()=ByteRGBValueC(255,0,0);
         }
      }
      Save("@X:left edge map overlaid on reference image", overlay, "", VERBOSE);
   }

   // Define region of interest in images
   //  by default this is not necessary
   //  however this becomes useful to eliminate e.g. a black border around the image
   IndexRange2dC activeWindow(org, end);
   if(activeWindow.IsEmpty()) {
      activeWindow = rgbImages[refCamID].Frame();
   }

   // compute allowed disparity values
   ImageC<RangeVectorC> VHrange;
   ImageC<PairC<IntT> > BGrange;

   ImageC<RealT> depthMap(rgbRefImage.Frame());
   if(!VH) {
      ComputeIntersectionRange(VHrange, BGrange, triMaps, refTriMapMask, depthImgs, cameras, activeWindow, refCamID, yPitchMin, yPitchMax, depthStep, maxSurfDist);

      RealT flow = GraphCut(depthMap, VHrange, BGrange, rgbImages, bgmeanImages, depthImgs, cameras, refCamID, oCamIDs, fgVars, bgVarImages, edgeMap, depthStep, tolVis, useMean, winsize, kaSmoothMin, kaSmoothMax, kdSmoothMin, kdSmoothMax, tData, sigmaFG, sigmaBG, kData, maxBackground, noBG, activeWindow, useBestScore, VERBOSE_FG, VERBOSE_BG, testImgPix);
      cerr << "max flow = " << flow << endl;

   } else {
      // find first entrance point
      for(Array2dIter3C<RealT, RealT, ByteT> it(depthMap, depthImgs[refCamID], refTriMap, activeWindow); it; it++)
      {
         if(it.Data3()!=0) {
            it.Data1() = it.Data2();
         }
      }
   }

  {
   RealT dmax = 0.0;
   for(Array2dIterC<RealT> it(depthMap); it; it++)
   {
      if(it.Data()>dmax)
         dmax = it.Data();
   }
   ImageC<ByteT> depthMap2(rgbRefImage.Frame());
   for(Array2dIter2C<ByteT, RealT> it(depthMap2, depthMap); it; it++)
   {
      it.Data1() = Round(it.Data2()/dmax*255);
   }
   // Save depth map
   if(!outputDEPTH.IsEmpty()) {
      sprintf(pFile,outputDEPTH.chars(),refCamID,frameNumber);
      if(!Save(pFile, depthMap2, "", VERBOSE)) 
      {
         cerr << "Failed to save output depthmap of foreground reconstruction to " << pFile << "." << endl;
         return 1;
      }
   }

   if(VERBOSE) {
      Save("@X:depth map", depthMap, "", VERBOSE);

      ImageC<ByteRGBValueC> segFGImg(rgbRefImage.Frame());
      segFGImg.Fill(BGPixel);
      for(Array2dIter3C<ByteRGBValueC, ByteRGBValueC, RealT> itimg(segFGImg, rgbRefImage, depthMap); itimg; itimg++) {
         if(itimg.Data3()>0.0) {
            itimg.Data1() = itimg.Data2();
         }
      }
      Save("@X:segmented foreground image", segFGImg, "", VERBOSE);
   }

   if(!outputMESH.IsEmpty()) {
      // Build mesh of foreground reconstruction
      BuildMesh(depthMap, triMaps, cameras, activeWindow, refCamID, depthStep, smooth, rgbRefImage, embedColour, frameNumber, outputMESH);
   }

  }

   if(REFINE_BOUNDARY) {
      // run refine matting on player boundary strips

      // build trimap from depthmap
      ImageC<ByteT> refinementTriMap(depthMap.Frame());
      refinementTriMap.Fill(0);
      ImageC<ByteT> tmpTriMap2(depthMap.Frame());
      tmpTriMap2.Fill(0);
      ImageC<ByteT> tmpTriMap3(depthMap.Frame());
      tmpTriMap3.Fill(0);
      for(Array2dIter3C<RealT, ByteT, ByteT> it(depthMap, tmpTriMap2, refinementTriMap); it; it++) {
         if(it.Data1()>0) {
            it.Data2() = 255;
            it.Data3() = 255;
         }
      }
      ImageC<ByteT> refinementKernel1(IndexRange2dC(-dilate_ref,dilate_ref,-dilate_ref,dilate_ref));
      refinementKernel1.Fill(127);
      BinaryDilate(tmpTriMap2, refinementKernel1, tmpTriMap3, ByteT(127));
      for(Array2dIter4C<ByteT, ByteT, ByteT, bool> it(tmpTriMap2, tmpTriMap3, refinementTriMap, refTriMapMask); it; it++) {
         if(it.Data1() == 0 && it.Data2()==127 && it.Data4()==true) {
            it.Data3() = 127;
         }
      }
      ImageC<ByteT> refinementKernel2(IndexRange2dC(-erode_ref,erode_ref,-erode_ref,erode_ref));
      refinementKernel2.Fill(127);
      tmpTriMap3.Fill(0);
      BinaryErode(tmpTriMap2, refinementKernel2, tmpTriMap3, ByteT(127));
      for(Array2dIter3C<ByteT, ByteT, ByteT> it(tmpTriMap2, tmpTriMap3, refinementTriMap); it; it++) {
         if(it.Data1() == 255 && it.Data2()!=127) {
            it.Data3() = 127;
         }
      }

      if(VERBOSE) {
        Save("@X:Eroded silhouette", tmpTriMap3);
      }

      tmpTriMap3.Fill(0);
      BinarySkeletonisation(tmpTriMap2, tmpTriMap3, ByteT(255));

      if(VERBOSE) {
        Save("@X:Skeleton", tmpTriMap3);
      }

      BinaryPruning(tmpTriMap3.Copy(), tmpTriMap3, ByteT(255), pruning_ref);
      for(Array2dIter3C<ByteT, ByteT, ByteT> it(tmpTriMap2, tmpTriMap3, refinementTriMap); it; it++) {
         if(it.Data1() == 255 && it.Data2()==255) {
            it.Data3() = 255;
         }
      }

      if(VERBOSE) {
        Save("@X:Pruned skeleton", tmpTriMap3);
      }

      if(VERBOSE) {
         Save("@X:refinementTriMap", refinementTriMap);
      }

      ImageC<ByteRGBValueC> fgResult; // foreground image computation not implemented yet
      Compute2DMatte(rgbRefImage, refinementTriMap);

      for(Array2dIterC<ByteT> it(refinementTriMap); it; it++) {
         if(it.Data()<cuttingVal) {
            it.Data() = 0;
         }
      }

      if(!outputMATTE.IsEmpty()) {
         sprintf(pFile,outputMATTE.chars(),refCamID,frameNumber);
         if(!Save(pFile, refinementTriMap, "", VERBOSE))
         {
            cerr << "Failed to save output alpha matte to " << pFile << "." << endl;
            return 1;
         }
      }
      if(VERBOSE) {
         Save("@X:Refined alpha matte", refinementTriMap);
      }

      if(!outputFG.IsEmpty()) {
         sprintf(pFile,outputFG.chars(),refCamID,frameNumber);
         if(!Save(pFile, fgResult, "", VERBOSE))
         {
            cerr << "Failed to save output foreground to " << pFile << "." << endl;
            return 1;
         }
      }
      if(VERBOSE) {
         Save("@X:Estimated foreground", fgResult);
      }

     if(!outputDEPTH2.IsEmpty() || !outputMESH2.IsEmpty()) {

      triMaps[refCamID] = refinementTriMap;

      ImageC<RealT> depthMap(rgbRefImage.Frame());
      if(!VH2) {
         ComputeIntersectionRange(VHrange, BGrange, triMaps, refTriMapMask, depthImgs, cameras, activeWindow, refCamID, yPitchMin, yPitchMax, depthStep, maxSurfDist);

         RealT flow = GraphCut(depthMap, VHrange, BGrange, rgbImages, bgmeanImages, depthImgs, cameras, refCamID, oCamIDs, fgVars, bgVarImages, edgeMap, depthStep, tolVis, useMean, winsize, kaSmoothMin, kaSmoothMax, kdSmoothMin, kdSmoothMax, tData, sigmaFG, sigmaBG, kData, maxBackground, true, activeWindow, useBestScore, VERBOSE_FG, VERBOSE_BG, testImgPix);
         cerr << "max flow = " << flow << endl;
      } else {
         // find first entrance point
         for(Array2dIter3C<RealT, RealT, ByteT> it(depthMap, depthImgs[refCamID], refinementTriMap, activeWindow); it; it++)
         {
            if(it.Data3()!=0) {
               it.Data1() = it.Data2();
            }
         }
      }

      RealT dmax = 0.0;
      for(Array2dIterC<RealT> it(depthMap); it; it++)
      {
         if(it.Data()>dmax)
            dmax = it.Data();
      }
      ImageC<ByteT> depthMap2(rgbRefImage.Frame());
      for(Array2dIter2C<ByteT, RealT> it(depthMap2, depthMap); it; it++)
      {
         it.Data1() = Round(it.Data2()/dmax*255);
      }
      // Save depth map
      if(!outputDEPTH2.IsEmpty()) {
         sprintf(pFile,outputDEPTH2.chars(),refCamID,frameNumber);
         if(!Save(pFile, depthMap2, "", VERBOSE)) 
         {
            cerr << "Failed to save output depthmap of foreground reconstruction to " << pFile << "." << endl;
            return 1;
         }
      }

      if(VERBOSE) {
         Save("@X:depth map 2", depthMap, "", VERBOSE);

         ImageC<ByteRGBValueC> segFGImg(rgbRefImage.Frame());
         segFGImg.Fill(BGPixel);
         for(Array2dIter3C<ByteRGBValueC, ByteRGBValueC, RealT> itimg(segFGImg, rgbRefImage, depthMap); itimg; itimg++) {
            if(itimg.Data3()>0.0) {
               itimg.Data1() = itimg.Data2();
            }
         }
         Save("@X:segmented foreground image 2", segFGImg, "", VERBOSE);
      }

      if(!outputMESH2.IsEmpty()) {
         // Build mesh of foreground reconstruction
         BuildMesh(depthMap, triMaps, cameras, activeWindow, refCamID, depthStep, smooth, rgbRefImage, embedColour, frameNumber, outputMESH2, refinementTriMap);
      }
     }
   }

   // Save stats
   if(!outputSTATS.IsEmpty()) {
      SArray1dC<ImageC<RealRGBValueC> > warpedImages(numCams);
      SArray1dC<ImageC<RealT> > errorMaps(numCams);
      ProjectToReference(rgbImages, depthMap, cameras, refCamID, oCamIDs, warpedImages, errorMaps, activeWindow, frameNumber, outputSTATS, useRobustStats, useIsotropic, VERBOSE);
   }

   return 0;
}


// Compute edge maps
bool
ComputeEdgeMap(const ImageC<ByteRGBValueC>& rgbImage,
               ImageC<bool>& edgeMap,
               RealT minHyst,
               RealT maxHyst,
               bool useSobel,
               RealT omegaDeriche,
               RealT alphaDeriche)
{

   ImageC<RealT> rdiff(rgbImage.Frame());
   ImageC<RealT> cdiff(rgbImage.Frame());
   for(Array2dIter3C<ByteRGBValueC, RealT, RealT> it(rgbImage, rdiff, cdiff); it; it++) {
      Index2dC rneighbour = it.Index()+Index2dC(1,0);
      if(rgbImage.Contains(rneighbour)) {
         Vector3dC v1(it.Data1().Red(),it.Data1().Green(),it.Data1().Blue());
         Vector3dC v2(rgbImage[rneighbour].Red(),rgbImage[rneighbour].Green(),rgbImage[rneighbour].Blue());
         it.Data2() = v1.EuclidDistance(v2);
      } else {
         it.Data2() = 0.0;
      }
      Index2dC cneighbour = it.Index()+Index2dC(0,1);
      if(rgbImage.Contains(cneighbour)) {
         Vector3dC v1(it.Data1().Red(),it.Data1().Green(),it.Data1().Blue());
         Vector3dC v2(rgbImage[cneighbour].Red(),rgbImage[cneighbour].Green(),rgbImage[cneighbour].Blue());
         it.Data3() = v1.EuclidDistance(v2);
      } else {
         it.Data3() = 0.0;
      }
   }

   // Convert image to grayscale
   ImageC<RealT> gsImage(rgbImage.Frame());
   for(Array2dIter2C<ByteRGBValueC, RealT> it(rgbImage, gsImage); it; it++) {
      it.Data2() = 0.3*it.Data1().Red() + 0.59*it.Data1().Green() + 0.11*it.Data1().Blue();
   }

   // Apply edge filter
   ImageC<RealT> edgeDr;
   ImageC<RealT> edgeDc;
   if(useSobel) {
      EdgeSobelC<RealT,RealT> edgeDet;
      edgeDet.Apply(gsImage,edgeDr,edgeDc);
   } else {
      EdgeDericheC edgeDet(omegaDeriche, alphaDeriche);
      edgeDet.Apply(gsImage,edgeDr,edgeDc);
   }

   // Compute edge magnitude
   ImageC<RealT> edgeMag;
   SqrCompositionC sqrCompose;
   sqrCompose.Apply(edgeDr,edgeDc,edgeMag);
   RavlAssert(edgeDr.Frame().Area() > 0);
   RavlAssert(edgeDc.Frame().Area() > 0);
   RavlAssert(edgeMag.Frame().Area() > 0);

   // Apply non maxima suppression
   ImageC<RealT> nonMax;
   RealT mean;
   IntT edgels;
   EdgeNonMaxSuppressionC nonMaxSup;
   nonMaxSup.Apply(edgeDr,edgeDc,edgeMag,nonMax,mean,edgels);

   // Apply hysteresis thresholding
   EdgeLinkC edgeLink = HysteresisThreshold(nonMax,minHyst,maxHyst);

   // Link edges
   SArray1dC<EdgelC> edges = edgeLink.ArrayOfEdgels(edgeDr,edgeDc,nonMax);

   // Save as an image
   edgeMap = ImageC<bool>(rgbImage.Frame());
   edgeMap.Fill(false);
   for(SArray1dIterC<EdgelC> it(edges);it;it++)
      edgeMap[it->At()] = true;

   return true;
}

// Determine the projection of the point Pt into the images with IDs CamIDs
bool
ComputeProjections(SArray1dC<Vector2dC>& imgPts,
                   const Vector3dC& Pt,
                   const PinholeCameraArrayC& cameras,
                   const SArray1dC<IntT>& camIDs)
{
   imgPts = SArray1dC<Vector2dC>(camIDs.Size());
   for(Array1dIter2C<IntT, Vector2dC> it(camIDs, imgPts); it; it++) {
      cameras[it.Data1()].Project(it.Data2(), Pt);
   }

   return true;
}

// Determine the projection of the point Pt into the images with IDs CamIDs,
//  partitioning cameras into foreground and background
bool
ComputeProjections(SArray1dC<Vector2dC>& imgPts,
                   DListC<IntT>& fgCamIDs,
                   DListC<IntT>& bgCamIDs,
                   DListC<IntT>& ocCamIDs,
                   const Vector3dC& Pt,
                   const PinholeCameraArrayC& cameras,
                   const SArray1dC<IntT>& camIDs,
                   const SArray1dC<ImageC<ByteT> >& triMaps,
                   const IndexRange2dC activeWindow)
{
   imgPts = SArray1dC<Vector2dC>(camIDs.Size());
   fgCamIDs.Empty(); bgCamIDs.Empty(); ocCamIDs.Empty();
   for(Array1dIter2C<IntT, Vector2dC> it(camIDs, imgPts); it; it++) {
      // compute projection
      cameras[it.Data1()].Project(it.Data2(), Pt);

      // check layer the projection belongs to
      Index2dC ind(Round(it.Data2().Y()),Round(it.Data2().X()));
      if(!activeWindow.Contains(ind)) {
         ocCamIDs+=it.Data1();
      } else if(triMaps[it.Data1()][ind]==0) {
         bgCamIDs+=it.Data1();
      } else {
         fgCamIDs+=it.Data1();
      }
   }

   return true;
}

// Determines if a 3D point is within the visual hull by looking at
//  its projection in the visual hull silhouettes
bool
IsWithinVisualHull(const Vector3dC& Pt,
                   const PinholeCameraArrayC& cameras,
                   const SArray1dC<ImageC<ByteT> >& triMaps,
                   const IndexRange2dC activeWindow,
                   bool verbose)
{
   SArray1dC<Vector2dC> imgPts;
   DListC<IntT> fgCamIDs, bgCamIDs, ocCamIDs;
   SArray1dC<IntT> camIDs(cameras.Size());
   for(Array1dIterC<IntT> it(camIDs); it; it++) {
      it.Data() = it.Index().V();
   }
   ComputeProjections(imgPts, fgCamIDs, bgCamIDs, ocCamIDs, Pt, cameras, camIDs, triMaps, activeWindow);

   if(fgCamIDs.IsEmpty() || !bgCamIDs.IsEmpty()) {
      return false;
   } else {
      return true;
   }
}

// Determine the intersection range for camera image
bool
ComputeIntersectionRange(ImageC<RangeVectorC>& VHrange,
                         ImageC<PairC<IntT> >& BGrange,
                         const SArray1dC<ImageC<ByteT> >& triMaps,
                         const ImageC<bool>& refTriMapMask,
                         const SArray1dC<ImageC<RealT> >& depthImgs,
                         const PinholeCameraArrayC& cameras,
                         const IndexRange2dC activeWindow,
                         IntT refCamID,
                         RealT yPitchMin,
                         RealT yPitchMax,
                         RealT step,
                         IntT maxSurfDist)
{

   const ImageC<ByteT>& refTriMap = triMaps[refCamID];
   const ImageC<RealT>& refDepthImg = depthImgs[refCamID];
   VHrange = ImageC<RangeVectorC>(refTriMap.Frame());
   BGrange = ImageC<PairC<IntT> >(refTriMap.Frame());

   // set invalid ranges for foreground and background
   PairC<IntT> invalidPair(invalidIndex, invalidIndex);
   BGrange.Fill(invalidPair);
   RangeVectorC invalidRangeVector; invalidRangeVector.Fill(invalidPair);
   VHrange.Fill(invalidRangeVector);

   // Get reference camera
   const PinholeCameraC& refCam = cameras[refCamID];
   Vector3dC refCamOrg;
   refCam.Origin(refCamOrg);

   // indicates if RangeVectorC is of sufficient size or not (true means some segments are missing)
   bool overflow = false;

   // Iterate over the reference camera trimap setting the index range for each pixel
   Array2dIter4C<ByteT, RangeVectorC, PairC<IntT>, bool> itimg(refTriMap, VHrange, BGrange, refTriMapMask, activeWindow);
   for (; itimg; itimg++)
   {
      ByteT& pixVal = itimg.Data1();
      RangeVectorC& VHrng = itimg.Data2();
      PairC<IntT>& BGrng = itimg.Data3();
      bool& pixMask = itimg.Data4();

      Vector2dC imgPt(itimg.Index().Col(), itimg.Index().Row());
      Index2dC imgPix(Round(imgPt.Y()),Round(imgPt.X()));

      // Determine the 3D direction of the ray
      Vector3dC dir;
      refCam.ProjectInverseDirection(dir, imgPt);
      dir.MakeUnit();

      // Determine depth of intersection with background (it is assumed the background is the plane Y=yPitchMax)
      RealT depthPitchMax = (yPitchMax-refCamOrg.Y())/dir.Y();

      // Get the depth of first intersection with visual hull from the depth buffer image
      RealT depthInter = refDepthImg[imgPix];

      // Set the range if this is not a background pixel,
      //  and there is a valid depth, i.e. the point is in front of the camera and above pitch level
      if (pixVal != 0 && depthInter > 0.0 && depthInter < depthPitchMax)
      {
         bool inside = false;
         IntT entranceDepth = 0, exitDepth, surfDist = 0;
         UIntT ind = 0;

         for(IntT depth=Ceil(depthInter/step); depth<Floor(depthPitchMax/step); depth++) {
            // Check if point belongs to visual hull
            Vector3dC pt0 = refCamOrg + dir*depth*step;
            bool isWithinVH = IsWithinVisualHull(pt0, cameras, triMaps, activeWindow);

            if(!inside && isWithinVH) {
               // entrance point found
               inside = true;
               entranceDepth = depth;
               surfDist = 0;
            } else if(inside && isWithinVH) {
               surfDist++;
               if(surfDist==maxSurfDist+1) {
                  exitDepth = depth-1;
                  // exit point found
                  VHrng[ind] = PairC<IntT>(entranceDepth, exitDepth);
                  ind++;
                  if(ind==VHrange.Size()) {
                     overflow = true;
                     break;
                  }
               }
            } else if(inside && !isWithinVH) {
               inside = false;
               if(surfDist<=maxSurfDist) {
                  exitDepth = depth-1;
                  // exit point found
                  VHrng[ind] = PairC<IntT>(entranceDepth, exitDepth);
                  ind++;
                  if(ind==VHrange.Size()) {
                     overflow = true;
                     break;
                  }
               }
            }

         }
         if(inside && (surfDist<=maxSurfDist)) {
            // close last segment at limit if needed
            inside = false;
            exitDepth = Floor(depthPitchMax/step)-1;
            VHrng[ind] = PairC<IntT>(entranceDepth, exitDepth);
         }
      }

      // Set the range if this is not a foreground pixel
      if (pixVal!=255 && pixMask==true)
      {
         // Determine exit point with background (it is assumed the outer background plane is the plane Y=YpitchMin)
         RealT depthPitchMin = (yPitchMin-refCamOrg.Y())/dir.Y();
         //BGrng = PairC<IntT>(Floor(depthPitchMax/step), Floor(depthPitchMin/step));
         BGrng = PairC<IntT>(99999999, 99999999);
      }
   }
   if(overflow) {
      cerr << "Warning: some intersections with visual hull may not have been considered. To guarantee all intestections are considered, increase maximum of segments allocated." << endl;
   }

   return true;

}

// Determine the visibility of the point Pt in the images with IDs CamIDs
bool
ComputeVisibility(SArray1dC<IntT>& visCamIDs,
                  const SArray1dC<ImageC<RealT> >& depthImgs,
                  const Vector3dC& Pt,
                  const PinholeCameraArrayC& cameras,
                  const SArray1dC<IntT>& camIDs,
                  RealT tolVis)
{
   SArray1dC<Vector2dC> imgPts;
   ComputeProjections(imgPts, Pt, cameras, camIDs);

   SArray1dC<bool> vis(camIDs.Size());
   IntT count = 0;
   for(Array1dIter3C<IntT, Vector2dC, bool> it(camIDs, imgPts, vis); it; it++) {
      // Compute current camera centre
      Vector3dC camOrg;
      cameras[it.Data1()].Origin(camOrg);

      // Compute distance from current camera centre
      RealT d = camOrg.EuclidDistance(Pt);

      // Compare against depth map
      Index2dC ImgPix(it.Data2().Y(), it.Data2().X());
      if(depthImgs[it.Data1()].Contains(ImgPix) &&
           (Abs(depthImgs[it.Data1()][ImgPix] - d)<tolVis || depthImgs[it.Data1()][ImgPix]<=0)) {
         it.Data3() = true;
         count++;
      } else {
         it.Data3() = false;
      }
   }

   visCamIDs = SArray1dC<IntT>(count);
   IntT ind = 0;
   for(Array1dIterC<bool> it(vis); it; it++) {
      if(it.Data()==true) {
         visCamIDs[ind] = camIDs[it.Index()];
         ind++;
      }
   }

   return true;
}

// Return smoothness cost adapted according to Euclidian distance
bool SmoothnessCost(RealT& score,
                    bool Edge1,
                    bool Edge2,
                    RealT kSmoothMin,
                    RealT kSmoothMax)
{
   score = (Edge1 || Edge2) ? kSmoothMin : kSmoothMax;

   return true;
}

// Return the sum of squarred difference between a reference image window and other image windows
//  using background information
// Returns false if the number of background pixels exceeds the tolerance
bool ComputeFGScoreMean(RealT& score,
                   const ImageC<ByteRGBValueC>& refPatch,
                   const SArray1dC<ImageC<ByteRGBValueC> >& oPatch,
                   const ImageC<ByteRGBValueC>& refBGPatch,
                   const SArray1dC<ImageC<ByteRGBValueC> >& oBGPatch,
                   RealT tData,
                   RealT kData,
                   SArray1dC<Vector3dC> fgVars,
                   Vector3dC refbgVar,
                   SArray1dC<Vector3dC> obgVars,
                   IntT maxBackground,
                   bool VERBOSE)
{
   IntT count = 0;
   score = 0.0;

   if(refPatch.IsEmpty()) {
      score = maxScore;
      return false;
   }

   RealT scorei;
   for(Array1dIter4C<ImageC<ByteRGBValueC>, ImageC<ByteRGBValueC>, Vector3dC, Vector3dC> it(oPatch, oBGPatch, fgVars, obgVars); it; it++) {
      if(!it.Data1().IsEmpty() && SSDCorrelation(scorei, refPatch, it.Data1(), refBGPatch, it.Data2(), tData, kData, it.Data3(), PairC<Vector3dC>(refbgVar,it.Data4()), maxBackground, VERBOSE))
      {
         if(VERBOSE) {
            cerr << "Partial score: " << scorei << endl;
         }
         score += scorei;
         count++;
      }
   }

   if(count!=0) {
      score/=count;
   } else {
      score = maxScore;
      return false;
   }

   if(VERBOSE) {
      cerr << "Mean score: " << score << endl;
   }

   return true;
}

// Return the sum of squarred difference between a reference image window and other image windows
//  using background information
// Returns false if the number of background pixels exceeds the tolerance
bool ComputeFGScoreBest(RealT& score,
                   const ImageC<ByteRGBValueC>& refPatch,
                   const SArray1dC<ImageC<ByteRGBValueC> >& oPatch,
                   const ImageC<ByteRGBValueC>& refBGPatch,
                   const SArray1dC<ImageC<ByteRGBValueC> >& oBGPatch,
                   RealT tData,
                   RealT kData,
                   SArray1dC<Vector3dC> fgVars,
                   Vector3dC refbgVar,
                   SArray1dC<Vector3dC> obgVars,
                   IntT maxBackground,
                   bool VERBOSE)
{
   IntT count = 0;
   score = maxScore;

   if(refPatch.IsEmpty()) {
      score = maxScore;
      return false;
   }

   RealT scorei=1.0;
   for(Array1dIter4C<ImageC<ByteRGBValueC>, ImageC<ByteRGBValueC>, Vector3dC, Vector3dC> it(oPatch, oBGPatch, fgVars, obgVars); it; it++) {
      if(!it.Data1().IsEmpty() && SSDCorrelation(scorei, refPatch, it.Data1(), refBGPatch, it.Data2(), tData, kData, it.Data3(), PairC<Vector3dC>(refbgVar,it.Data4()), maxBackground, VERBOSE))
      {
         if(VERBOSE) {
            cerr << "Partial score: " << scorei << endl;
         }
         score = Min(scorei,score);
         count++;
      }
   }

   if(count==0) {
      return false;
   }

   if(VERBOSE) {
      cerr << "Best score: " << score << endl;
   }

   return true;
}


// Return the sum of squarred difference between two image windows
//  using background information
// Returns false if the number of background pixels exceeds the tolerance
bool SSDCorrelation(RealT& score,
                    const ImageC<ByteRGBValueC>& keyPatch,
                    const ImageC<ByteRGBValueC>& offsetPatch,
                    const ImageC<ByteRGBValueC>& keyBGPatch,
                    const ImageC<ByteRGBValueC>& offsetBGPatch,
                    RealT tData,
                    RealT kData,
                    Vector3dC fgVar,
                    PairC<Vector3dC> bgVars,
                    IntT maxBackground,
                    bool VERBOSE)
{
   // By default the score is high (bad match)
   score=maxScore;

   RealT FGscore, keyBGscore, offBGscore;
   if(SSDCorrelation(FGscore, keyPatch, offsetPatch, fgVar, maxBackground, false))
   {
      ImageC<ByteRGBValueC> keyPatch_center(keyPatch, IndexRange2dC(keyPatch.Frame().Center(), 0));
      SSDCorrelation(keyBGscore, keyPatch_center, keyBGPatch, bgVars.Data1(), maxBackground, false);
      ImageC<ByteRGBValueC> offsetPatch_center(offsetPatch, IndexRange2dC(offsetPatch.Frame().Center(), 0));
      SSDCorrelation(offBGscore, offsetPatch_center, offsetBGPatch, bgVars.Data2(), maxBackground, false);
      //score = (keyBGscore>tData && offBGscore>tData) ? FGscore : kData*FGscore;
      score = FGscore;
      score = (keyBGscore>tData) ? score : kData*score;
      score = (offBGscore>tData) ? score : kData*score;
   } else {
      return false;
   }

   if(VERBOSE) {
      cerr << "FG: " << FGscore << "(variance: " << fgVar << ")" << "  ";
      cerr << "BG1: " << keyBGscore << "(threshold: " << tData << ", variance: " << bgVars.Data1() << ")  ";
      cerr << "BG2: " << offBGscore << "(threshold: " << tData << ", variance: " << bgVars.Data2() << ")  ";
      cerr << "score: " << score << endl;
   }

   return true;
}

// Return the sum of squarred difference between two image windows
// Returns false if the number of background pixels exceeds the tolerance
bool SSDCorrelation(RealT& score,
                    const ImageC<ByteRGBValueC>& keyPatch,
                    const ImageC<ByteRGBValueC>& offsetPatch,
                    Vector3dC var,
                    IntT maxBackground,
                    bool VERBOSE)
{
   RavlAssertMsg(!keyPatch.IsEmpty(), "keyPatch cannot be empty.");
   RavlAssertMsg(!offsetPatch.IsEmpty(), "offsetPatch cannot be empty.");
   RavlAssertMsg(keyPatch.Size() == offsetPatch.Size(), "keyPatch and offsetPatch must be the same size.");

   // By default the score is high (bad match)
   score=maxScore;

   // Compute the statistics in the FG patches
   RealT sum = 0.0;
   IntT count = 0, back = 0;
   Array2dIter2C<ByteRGBValueC,ByteRGBValueC> it(keyPatch,offsetPatch);
   for (; it; it++) 
   {
      ByteRGBValueC kval = it.Data1();
      ByteRGBValueC oval = it.Data2();
      if (kval!=BGPixel && oval!=BGPixel)
      {
         count++;
         sum += Sqr(oval.Red()-kval.Red())/var.X()+Sqr(oval.Green()-kval.Green())/var.Y()+Sqr(oval.Blue()-kval.Blue())/var.Z();
      }
      else
      {
         back++;
         if (back > maxBackground) return false;
      }
   }

   if(count==0) {
      return false;
   }

   // return mean normalised between 0+eps and 1+eps
   score = sum/(count*3*Sqr(255));

   // clamp between 0 and maxScore
   score = Min(score, maxScore);

   if(VERBOSE) {
      cerr << "score: " << score << endl;
   }

   return true;
}

// Perform a graph-cut to maximise photoconsistency between views
//  imposing smoothness constaint based on edgeness
RealT
GraphCut(ImageC<RealT>& depthMap,
         const ImageC<RangeVectorC>& VHrange,
         const ImageC<PairC<IntT> >& BGrange,
         const SArray1dC<ImageC<ByteRGBValueC> >& rgbImages,
         const SArray1dC<ImageC<ByteRGBValueC> >& bgmeanImages,
         const SArray1dC<ImageC<RealT> >& depthImgs,
         const PinholeCameraArrayC& cameras,
         IntT refCamID,
         const SArray1dC<IntT>& oCamIDs,
         const SArray1dC<Vector3dC>& fgVars,
         const SArray1dC<ImageC<RealRGBValueC> >& bgVarImages,
         const ImageC<bool>& edgeMap,
         RealT step,
         RealT tolVis,
         bool useMean,
         IntT winsize,
         RealT kaSmoothMin,
         RealT kaSmoothMax,
         RealT kdSmoothMin,
         RealT kdSmoothMax,
         RealT tData,
         RealT sigmaFG,
         RealT sigmaBG,
         RealT kData,
         IntT maxBackground,
         bool noBG,
         IndexRange2dC activeWindow,
         bool useBestScore,
         bool FG_VERBOSE,
         bool BG_VERBOSE,
         Index2dC testImgPix)
{

   // We assume here that the capacities are normally in the interval [0,1]
   // so that InfiniteCap >> normal capacities
   Graph::captype InfiniteCap = 1e9;

   // set region of interest to the whole image if unspecified
   if(activeWindow.IsEmpty()) {
      activeWindow = rgbImages[refCamID].Frame();
   }

   IntT BGlayer;
   if(noBG) {
      BGlayer = 0;
   } else {
      BGlayer = 1;
   }

   // Create the graph structure
   Graph *g = new Graph();

   // Create a node for every point in the disparity space
   // Nodes are indexed as Index3dC(leftrow, leftcol, depthIndex)
   HashC< Index3dC, Graph::node_id > nVHlookup; // foreground node
   HashC< Index3dC, Graph::node_id > nBGlookup; // background node
   HashC< Index3dC, Graph::node_id > nEXlookup; // auxiliary node
   Array2dIter2C<RangeVectorC, PairC<IntT> > itl(VHrange, BGrange);
   IndexRange2dIterC itlpix(VHrange.Frame());
   RealT wka1, wka2;
   for (; itl; itl++, itlpix++)
   {
      // initialise with invalid nodes
      IntT firstCol = invalidIndex, lastCol = invalidIndex;

      // compute smoothness terms
      Index2dC n1(itlpix.Data().Row().V()-1, itlpix.Data().Col().V());
      if(rgbImages[refCamID].Contains(n1)) {
         ByteRGBValueC p0 = rgbImages[refCamID][itlpix.Data()];
         ByteRGBValueC p1 = rgbImages[refCamID][n1];
         Vector3dC v0(p0.Red(),p0.Green(),p0.Blue());
         Vector3dC v1(p1.Red(),p1.Green(),p1.Blue());
         wka1 = kaSmoothMin*Exp(-Sqr(v0.EuclidDistance(v1))/(2*Sqr(kaSmoothMax)));
      } else {
         wka1 = 0.0;
      }
      Index2dC n2(itlpix.Data().Row().V(), itlpix.Data().Col().V()-1);
      if(rgbImages[refCamID].Contains(n2)) {
         ByteRGBValueC p0 = rgbImages[refCamID][itlpix.Data()];
         ByteRGBValueC p2 = rgbImages[refCamID][n2];
         Vector3dC v0(p0.Red(),p0.Green(),p0.Blue());
         Vector3dC v2(p2.Red(),p2.Green(),p2.Blue());
         wka2 = kaSmoothMin*Exp(-Sqr(v0.EuclidDistance(v2))/(2*Sqr(kaSmoothMax)));
      } else {
         wka2 = 0.0;
      }
      // iterate through segments
      for(UIntT i=0; i<itl.Data1().Size()+BGlayer; i++) {
         // create reference to current look-up table
         HashC< Index3dC, Graph::node_id >& nlookup = (i<itl.Data1().Size()) ? nVHlookup : nBGlookup;
         // create reference to current segment of current layer
         PairC<IntT>& seg = (i<itl.Data1().Size()) ? itl.Data1()[i] : itl.Data2();

         // skip segment of empty
         if (seg.Data1() == invalidIndex || seg.Data2() == invalidIndex){
            continue;
         }

         // add previous node if missing (otherwise the first edge will not be assigned later)
         Index3dC v(itlpix.Data().Row().V(), itlpix.Data().Col().V(), seg.Data1()-1);
         Graph::node_id n;
         // if the previous node exists, it is the last node of previous segment and can only be part of FG
         if (!nVHlookup.Lookup(v,n)) {
            n = g->add_node();
            nEXlookup[ v ] = n;
         }

         // compute smoothness terms
         //RealT wka1 = kaSmoothMax;
         RealT wkd1 = kdSmoothMax;
         Index2dC n1(itlpix.Data().Row().V()-1, itlpix.Data().Col().V());
         if(edgeMap.Contains(n1)) {
            //SmoothnessCost(wka1, edgeMap[itlpix.Data()], edgeMap[n1], kaSmoothMin, kaSmoothMax);
            SmoothnessCost(wkd1, edgeMap[itlpix.Data()], edgeMap[n1], kdSmoothMin, kdSmoothMax);
         }

         //RealT wka2 = kaSmoothMax;
         RealT wkd2 = kdSmoothMax;
         Index2dC n2(itlpix.Data().Row().V(), itlpix.Data().Col().V()-1);
         if(edgeMap.Contains(n2)) {
            //SmoothnessCost(wka2, edgeMap[itlpix.Data()], edgeMap[n2], kaSmoothMin, kaSmoothMax);
            SmoothnessCost(wkd2, edgeMap[itlpix.Data()], edgeMap[n2], kdSmoothMin, kdSmoothMax);
         }
         // connect previous node to last node of previous segment with infinity capacity
         //  (otherwise there will be a break between source and sink)
         if (lastCol!=invalidIndex)
         {
            Index3dC vnn(itlpix.Data().Row().V(), itlpix.Data().Col().V(), lastCol);
            Graph::node_id nn;
            // the last node of previous segment can only be part of FG
            if (nVHlookup.Lookup(vnn, nn))
            {
               g -> add_edge(n, nn, InfiniteCap, InfiniteCap);
            }
         }

         if(i==itl.Data1().Size()) {
            // add penalty for changing layer
            IntT z1=0;
            Index2dC n1(itlpix.Data().Row().V()-1, itlpix.Data().Col().V());
            if(BGrange.Contains(n1)) {
               if(BGrange[n1].Data1()==invalidIndex) {
                  // find last exit point
                  for(UIntT k=0; k<VHrange[n1].Size(); k++) {
                     if(VHrange[n1][k].Data2() != invalidIndex) {
                        z1 = VHrange[n1][k].Data2();
                     }
                  }
               } else {
                  // find first entrance point
                  z1 = BGrange[n1].Data1()-1;
               }
               Index3dC vnn(n1.Row().V(), n1.Col().V(), z1);
               Graph::node_id nn;
               if (nVHlookup.Lookup(vnn, nn))
               {
                  g -> add_edge(n, nn, wka1, wka1);
               }
               else if (nBGlookup.Lookup(vnn, nn))
               {
                  g -> add_edge(n, nn, wka1, wka1);
               }
               else if (nEXlookup.Lookup(vnn, nn))
               {
                  g -> add_edge(n, nn, wka1, wka1);
               }
            }

            // add penalty for changing layer
            IntT z2=0;
            Index2dC n2(itlpix.Data().Row().V(), itlpix.Data().Col().V()-1);
            if(BGrange.Contains(n2)) {
               if(BGrange[n2].Data1()==invalidIndex) {
                  // find last exit point
                  for(UIntT k=0; k<VHrange[n2].Size(); k++) {
                     if(VHrange[n2][k].Data2() != invalidIndex) {
                        z2 = VHrange[n2][k].Data2();
                     }
                  }
               } else {
                  // find first entrance point
                  z2 = BGrange[n2].Data1()-1;
               }
               Index3dC vnn(n2.Row().V(), n2.Col().V(), z2);
               Graph::node_id nn;
               if (nVHlookup.Lookup(vnn, nn))
               {
                  g -> add_edge(n, nn, wka2, wka2);
               }
               else if (nBGlookup.Lookup(vnn, nn))
               {
                  g -> add_edge(n, nn, wka2, wka2);
               }
               else if (nEXlookup.Lookup(vnn, nn))
               {
                  g -> add_edge(n, nn, wka2, wka2);
               }
            }
         }

        // add each node of segment and create data and smoothness edges
         for (IntT depth = seg.Data1(); depth <= seg.Data2(); depth++)
         {
            // add current node
            Index3dC v(itlpix.Data().Row().V(), itlpix.Data().Col().V(), depth);
            Graph::node_id n;
            if (!nVHlookup.Lookup(v, n) && !nBGlookup.Lookup(v, n) && !nEXlookup.Lookup(v, n))
            {
               n = g->add_node();
               nlookup[ v ] = n;
            }

            // Add first neighbourhood edge to current
            Index3dC vnn(itlpix.Data().Row().V()-1, itlpix.Data().Col().V(), depth);
            Graph::node_id nn;
            if (nVHlookup.Lookup(vnn, nn))
            {
               g -> add_edge(n, nn, wkd1, wkd1);
            }
            else if (nBGlookup.Lookup(vnn, nn))
            {
               g -> add_edge(n, nn, wkd1, wkd1);
            }
            else if (nEXlookup.Lookup(vnn, nn))
            {
               g -> add_edge(n, nn, wkd1, wkd1);
            }

            // Add second neighbourhood edge to current
            vnn = Index3dC(itlpix.Data().Row().V(), itlpix.Data().Col().V()-1, depth);
            if (nVHlookup.Lookup(vnn, nn))
            {
               g -> add_edge(n, nn, wkd2, wkd2);
            }
            else if (nBGlookup.Lookup(vnn, nn))
            {
               g -> add_edge(n, nn, wkd2, wkd2);
            }
            else if (nEXlookup.Lookup(vnn, nn))
            {
               g -> add_edge(n, nn, wkd2, wkd2);
            }
         }
         if(firstCol==invalidIndex) {
            firstCol = seg.Data1()-1;
         }
         lastCol = seg.Data2();
      }
      if (firstCol != invalidIndex && lastCol != invalidIndex){
         Graph::node_id n;

         // Create source connection
         Index3dC s(itlpix.Data().Row().V(), itlpix.Data().Col().V(), firstCol);
         if (nVHlookup.Lookup(s, n))
         {
            g -> set_tweights(n, InfiniteCap, 0);
         }
         else if (nEXlookup.Lookup(s, n))
         {
            g -> set_tweights(n, InfiniteCap, 0);
         }

         // Create sink connection
         Index3dC t(itlpix.Data().Row().V(), itlpix.Data().Col().V(), lastCol);
         if (nVHlookup.Lookup(t, n))
         {
            g -> set_tweights(n, 0, InfiniteCap);
         }
         else if (nBGlookup.Lookup(t, n))
         {
            g -> set_tweights(n, 0, InfiniteCap);
         }
      }
   }

   // Calculate the photoconsistency score at every foreground node and place on the data edges
   // The cost is placed on the source side so that a cut edge corresponds to the node on the sink side
   HashIterC< Index3dC, Graph::node_id > itnVH(nVHlookup);
   for (; itnVH; itnVH++)
   {
      const Index3dC& v = itnVH.Key();
      Graph::node_id& n = itnVH.Data();
      Index2dC refImgPix(v.I(), v.J());

      // Create reference patch
      ImageC<ByteRGBValueC> refPatch,  refBGPatch;
      if (activeWindow.Contains(IndexRange2dC(refImgPix, winsize)))
      {
         refPatch = ImageC<ByteRGBValueC>(rgbImages[refCamID], IndexRange2dC(refImgPix, winsize));
      }
      else
      {
         refPatch = ImageC<ByteRGBValueC>(IndexRange2dC(refImgPix, winsize));
         refPatch.Fill(BGPixel);
         Array2dIter2C<ByteRGBValueC, ByteRGBValueC> itk(rgbImages[refCamID], refPatch, refPatch.Frame().ClipBy(activeWindow));
         for (; itk; itk++)
         {
            itk.Data2() = itk.Data1();
         }
      }

      if (activeWindow.Contains(refImgPix))
      {
         refBGPatch = ImageC<ByteRGBValueC>(bgmeanImages[refCamID], IndexRange2dC(refImgPix, 0));
      }

      // Get variance for reference image
      Vector3dC refbgVar(bgVarImages[refCamID][refImgPix].Red(),
                         bgVarImages[refCamID][refImgPix].Green(),
                         bgVarImages[refCamID][refImgPix].Blue());

      // Determine the 3D direction of the ray
      Vector3dC dir;
      Vector2dC refImgPt(refImgPix.Col(), refImgPix.Row());
      cameras[refCamID].ProjectInverseDirection(dir, refImgPt);
      dir.MakeUnit();

      // Determine 3D point P
      Vector3dC refCamOrg;
      cameras[refCamID].Origin(refCamOrg);
      Vector3dC Pt = refCamOrg + dir*v.K()*step;
      SArray1dC<Vector2dC> imgPts;
      ComputeProjections(imgPts, Pt, cameras, oCamIDs);

      SArray1dC<IntT> visCamIDs;
      // visibility handling
      if(depthImgs.Size()==0) {
         // do not worry about visibility
         visCamIDs = oCamIDs;
      } else {
         // Find closest entrance point to P on camera side
         IntT k0=v.K().V();
         Graph::node_id n0;
         while(nVHlookup.Lookup(Index3dC(v.I(), v.J(), IndexC(k0)), n0)) {
            k0--;
         }
         k0++;
         Vector3dC Pt0 = refCamOrg + dir*k0*step;

         // Calculate visibility
         ComputeVisibility(visCamIDs, depthImgs, Pt0, cameras, oCamIDs, tolVis);

         // if invisible
         if(visCamIDs.Size()==0) {
#if FALSE
            // Create data edges with maximum score
            Graph::node_id nn;
            if (nVHlookup.Lookup(Index3dC(v.I(), v.J(), v.K()-1), nn))
            {
               g -> add_edge(nn, n, maxScore, InfiniteCap);
            }
            else if (nBGlookup.Lookup(Index3dC(v.I(), v.J(), v.K()-1), nn))
            {
               g -> add_edge(nn, n, maxScore, InfiniteCap);
            }
            else if (nEXlookup.Lookup(Index3dC(v.I(), v.J(), v.K()-1), nn))
            {
               g -> add_edge(nn, n, maxScore, InfiniteCap);
            }

            continue;
#else
            // or use all cameras
            visCamIDs = oCamIDs;
#endif
         }
      }

      // Create other patches
      SArray1dC<ImageC<ByteRGBValueC> > visPatch(visCamIDs.Size());
      SArray1dC<ImageC<ByteRGBValueC> > visBGPatch(visCamIDs.Size());
      SArray1dC<Vector3dC> obgVars(visCamIDs.Size());
      for(Array1dIterC<IntT> it(visCamIDs); it; it++) {
         Index2dC ImgPix(imgPts[it.Index()].Y(), imgPts[it.Index()].X());
         if(!activeWindow.Contains(ImgPix)) {
            continue;
         }

         // Get variance for background image
         obgVars[it.Index()] = Vector3dC(bgVarImages[it.Data()][ImgPix].Red(),
                                         bgVarImages[it.Data()][ImgPix].Green(),
                                         bgVarImages[it.Data()][ImgPix].Blue());

         if (activeWindow.Contains(IndexRange2dC(ImgPix, winsize)))
         {
            visPatch[it.Index()] = ImageC<ByteRGBValueC>(rgbImages[it.Data()], IndexRange2dC(ImgPix, winsize));
         }
         else
         {
            visPatch[it.Index()] = ImageC<ByteRGBValueC>(IndexRange2dC(ImgPix, winsize));
            visPatch[it.Index()].Fill(BGPixel);
            Array2dIter2C<ByteRGBValueC, ByteRGBValueC> itk(rgbImages[it.Data()], visPatch[it.Index()], visPatch[it.Index()].Frame().ClipBy(rgbImages[it.Data()].Frame()));
            for (; itk; itk++)
            {
              itk.Data2() = itk.Data1();
            }
         }

         if (activeWindow.Contains(ImgPix))
         {
            visBGPatch[it.Index()] = ImageC<ByteRGBValueC>(bgmeanImages[it.Data()], IndexRange2dC(ImgPix, 0));
         }
      }

      // replace reference patch by mean (this is a symmetric matching score)
      if(useMean) {
         // compute mean patch
         SArray1dC<ImageC<ByteRGBValueC> > arefPatch(1);
         arefPatch[0] = refPatch.Copy();
         visPatch = arefPatch.Join(visPatch);
         ImageC<RealRGBValueC> meanPatch(refPatch.Frame());
         meanPatch.Fill(BGPixel);
         ImageC<IntT> count(refPatch.Frame());
         count.Fill(0);
         for(Array1dIterC<ImageC<ByteRGBValueC> > it(visPatch); it; it++) {
            for(Array2dIter3C<ByteRGBValueC, RealRGBValueC, IntT> itim(it.Data(), meanPatch, count); itim; itim++) {
               if(itim.Data1()!=BGPixel) {
                  itim.Data2().Red() += itim.Data1().Red();
                  itim.Data2().Green() += itim.Data1().Green();
                  itim.Data2().Blue() += itim.Data1().Blue();
                  itim.Data3()++;
               }
            }
         }
         refPatch = ImageC<ByteRGBValueC>(refPatch.Frame());
         for(Array2dIter3C<ByteRGBValueC, RealRGBValueC, IntT> itim(refPatch, meanPatch, count); itim; itim++) {
            if(itim.Data3()!=0) {
               itim.Data1().Red() = (IntT) (itim.Data2().Red()/itim.Data3());
               itim.Data1().Green() = (IntT) (itim.Data2().Green()/itim.Data3());
               itim.Data1().Blue() = (IntT) (itim.Data2().Blue()/itim.Data3());
            }
         }
         SArray1dC<ImageC<ByteRGBValueC> > arefBGPatch(1);
         arefBGPatch[0] = refBGPatch.Copy();
         visBGPatch = arefBGPatch.Join(visBGPatch);
         refBGPatch = ImageC<ByteRGBValueC>(refBGPatch.Frame());
         refBGPatch.Fill(BGPixel);
      }

      // Compute the correlation score
      RealT score;
      if(useBestScore) {
         ComputeFGScoreBest(score, refPatch, visPatch, refBGPatch, visBGPatch, tData, kData, fgVars, refbgVar, obgVars,  maxBackground, FG_VERBOSE || refImgPix==testImgPix);
      } else {
         ComputeFGScoreMean(score, refPatch, visPatch, refBGPatch, visBGPatch, tData, kData, fgVars, refbgVar, obgVars, maxBackground, FG_VERBOSE || refImgPix==testImgPix);
      }
      score = score/(sigmaFG*sigmaFG);

      // Create data edges
      Graph::node_id nn;
      if (nVHlookup.Lookup(Index3dC(v.I(), v.J(), v.K()-1), nn))
      {
         g -> add_edge(nn, n, score, InfiniteCap);
      }
      else if (nBGlookup.Lookup(Index3dC(v.I(), v.J(), v.K()-1), nn))
      {
         g -> add_edge(nn, n, score, InfiniteCap);
      }
      else if (nEXlookup.Lookup(Index3dC(v.I(), v.J(), v.K()-1), nn))
      {
         g -> add_edge(nn, n, score, InfiniteCap);
      }

      if(FG_VERBOSE || refImgPix==testImgPix) {
         ImageC<ByteRGBValueC> refOverlay(rgbImages[refCamID].Copy());
         ImageC<ByteRGBValueC> border(2*(2*winsize+3),(1+visPatch.Size())*(2*winsize+3)); border.Fill(ByteRGBValueC(255,0,0));
         refOverlay.SetSubArray(Index2dC(0,0), border);
         refOverlay.SetSubArray(Index2dC(winsize+1,winsize+1)-refImgPix, refPatch);
         refOverlay.SetSubArray(Index2dC(1*(2*winsize+3)+winsize+1,winsize+1)-refImgPix, refBGPatch);
         for(UIntT i=0; i<visPatch.Size(); i++) {
            Index2dC imgPix = (visPatch[i].Frame().TopLeft()+visPatch[i].Frame().BottomRight())/2;
            refOverlay.SetSubArray(Index2dC(winsize+1,(i+1)*(2*winsize+3)+winsize+1)-imgPix, visPatch[i]);
            refOverlay.SetSubArray(Index2dC((2*winsize+3)+winsize+1,(i+1)*(2*winsize+3)+winsize+1)-imgPix, visBGPatch[i]);
         }
         DrawCross(refOverlay, ByteRGBValueC(255,0,0), refImgPix, 3);
         Save("@X:refOverlay", refOverlay, "", true);

         for(UIntT j=0; j<visPatch.Size(); j++) {
            ImageC<ByteRGBValueC> offOverlay(rgbImages[visCamIDs[j]].Copy());
            ImageC<ByteRGBValueC> border(2*(2*winsize+3),(1+visPatch.Size())*(2*winsize+3)); border.Fill(ByteRGBValueC(255,0,0));
            offOverlay.SetSubArray(Index2dC(0,0), border);
            offOverlay.SetSubArray(Index2dC(winsize+1,winsize+1)-refImgPix, refPatch);
            offOverlay.SetSubArray(Index2dC(1*(2*winsize+3)+winsize+1,winsize+1)-refImgPix, refBGPatch);
            for(UIntT i=0; i<visPatch.Size(); i++) {
               Index2dC imgPix = (visPatch[i].Frame().TopLeft()+visPatch[i].Frame().BottomRight())/2;
               offOverlay.SetSubArray(Index2dC(winsize+1,(i+1)*(2*winsize+3)+winsize+1)-imgPix, visPatch[i]);
               offOverlay.SetSubArray(Index2dC((2*winsize+3)+winsize+1,(i+1)*(2*winsize+3)+winsize+1)-imgPix, visBGPatch[i]);
            }
            Index2dC offImgPix(imgPts[j].Y(), imgPts[j].X());
            DrawCross(offOverlay, ByteRGBValueC(255,0,0), offImgPix, 3);
            StringC stringOverlay = "@X:offOverlay"+StringC(j);
            Save(stringOverlay, offOverlay, "", true);
         }

         char tmp;
         cin >> tmp;
      }
   }

   // Calculate the photoconsistency score at every background node and place on the data edges
   // The cost is placed on the source side so that a cut edge corresponds to the node on the sink side
   HashIterC< Index3dC, Graph::node_id > itnBG(nBGlookup);
   for (; itnBG; itnBG++)
   {
      const Index3dC& v = itnBG.Key();
      Graph::node_id& n = itnBG.Data();
      Index2dC refImgPix(v.I(), v.J());
      // Create reference patch
      ImageC<ByteRGBValueC> refPatch,  refBGPatch;
      if (activeWindow.Contains(refImgPix))
      {
         refPatch = ImageC<ByteRGBValueC>(rgbImages[refCamID], IndexRange2dC(refImgPix, 0));
         refBGPatch = ImageC<ByteRGBValueC>(bgmeanImages[refCamID], IndexRange2dC(refImgPix, 0));
      }
      // Get variance for reference image
      Vector3dC refbgVar(bgVarImages[refCamID][refImgPix].Red(),
                         bgVarImages[refCamID][refImgPix].Green(),
                         bgVarImages[refCamID][refImgPix].Blue());

      // Compute the correlation score
      RealT score;
      SSDCorrelation(score, refPatch, refBGPatch, refbgVar, maxBackground, BG_VERBOSE || refImgPix==testImgPix);
      score = score/(sigmaBG*sigmaBG);

      // Create data edges
      Graph::node_id nn;
      if (nVHlookup.Lookup(Index3dC(v.I(), v.J(), v.K()-1), nn))
      {
         g -> add_edge(nn, n, score, InfiniteCap);
      }
      else if (nBGlookup.Lookup(Index3dC(v.I(), v.J(), v.K()-1), nn))
      {
         g -> add_edge(nn, n, score, InfiniteCap);
      }
      else if (nEXlookup.Lookup(Index3dC(v.I(), v.J(), v.K()-1), nn))
      {
         g -> add_edge(nn, n, score, InfiniteCap);
      }

      if(BG_VERBOSE || refImgPix==testImgPix) {
         ImageC<ByteRGBValueC> refOverlay(rgbImages[refCamID].Copy());
         ImageC<ByteRGBValueC> border(6,3); border.Fill(ByteRGBValueC(255,0,0));
         refOverlay.SetSubArray(Index2dC(0,0), border);
         refOverlay.SetSubArray(Index2dC(1,1)-refImgPix, refPatch);
         refOverlay.SetSubArray(Index2dC(4,1)-refImgPix, refBGPatch);
         DrawCross(refOverlay, ByteRGBValueC(255,0,0), refImgPix, 3);
         Save("@X:refOverlay", refOverlay, "", true);

         char tmp;
         cin >> tmp;
      }
   }

   // Perform max-flow
   Graph::flowtype flow = g->maxflow();

   // Extract correspondence as the sink side nodes of the cut
   depthMap = ImageC<RealT>(VHrange.Frame());
   Array2dIterC<RealT> itlm(depthMap);
   for (itl.First(), itlpix.First(); itl; itl++, itlpix++, itlm++)
   {
      itlm.Data() = 0;
      for(UIntT i=0; i<itl.Data1().Size(); i++) {
         HashC< Index3dC, Graph::node_id >& nlookup = (i<itl.Data1().Size()) ? nVHlookup : nBGlookup;
         PairC<IntT>& seg = (i<itl.Data1().Size()) ? itl.Data1()[i] : itl.Data2();
         if (seg.Data1() == invalidIndex || seg.Data2() == invalidIndex){
            continue;
         }
         for (IntT depth = seg.Data1(); depth <= seg.Data2(); depth++)
         {
            Graph::node_id n = nlookup[ Index3dC(itlpix.Data().Row().V(), itlpix.Data().Col().V(), depth-1) ];
            if (!n)
            {
               n = nEXlookup[ Index3dC(itlpix.Data().Row().V(), itlpix.Data().Col().V(), depth-1) ];
            }
            if (!n && i==itl.Data1().Size())
            {
               n = nVHlookup[ Index3dC(itlpix.Data().Row().V(), itlpix.Data().Col().V(), depth-1) ];
            }
            Graph::node_id nn = nlookup[ Index3dC(itlpix.Data().Row().V(), itlpix.Data().Col().V(), depth) ];
            if (!n || !nn)
            {
               cerr << "WARNING: some foreground nodes are missing!" << endl;
            }
            if (n && nn && g->what_segment(n) == Graph::SOURCE && g->what_segment(nn) == Graph::SINK)
            {
               itlm.Data()=depth*step;
               break;
            }
         }
      }
   }

   delete g;

   return (RealT)flow;

}


// Builds a mesh of the reconstruction
//  only points located at a minimum distance depthLimit from the camera are included in the mesh
bool
BuildMesh(const ImageC<RealT>& depthMap,
          const SArray1dC<ImageC<ByteT> >& triMaps,
          const PinholeCameraArrayC& cameras,
          const IndexRange2dC activeWindow,
          IntT refCamID,
          RealT step,
          IntT smooth,
          const ImageC<ByteRGBValueC>& fgRefImage,
          bool embedColour,
          IntT frameNumber,
          StringC output,
          const ImageC<ByteT>& alphaMatte)
{

   // Get reference camera
   const PinholeCameraC& refCam = cameras[refCamID];

   // Get camera origin
   Vector3dC camOrg;
   refCam.Origin(camOrg);

   // Process each pixel
   UIntT nv = 0;
   Array1dC<Vector3dC> v(depthMap.Size());
   ImageC<UIntT> vlookup(depthMap.Frame());
   vlookup.Fill(0);
   for(Array2dIter2C<RealT, UIntT> it(depthMap, vlookup); it; it++) {

      if(it.Data1()>0) {
         // Get image point coordinates
         Vector2dC imgPt(it.Index().Col(), it.Index().Row());

         // Determine the 3D direction of the ray
         Vector3dC dir;
         refCam.ProjectInverseDirection(dir, imgPt);
         dir.MakeUnit();

         // Compute 3D point coordinates
         Vector3dC pt = camOrg + dir*it.Data1();

         // Store data
         v[nv] = pt;
         it.Data2() = nv;
         nv++;
      }
   }

   // COMPUTE CONNECTIVITY
   ImageC<bool> rnimg(depthMap.Frame()); rnimg.Fill(false);
   ImageC<bool> bnimg(depthMap.Frame()); bnimg.Fill(false);
   ImageC<bool> brnimg(depthMap.Frame()); brnimg.Fill(false);
   ImageC<bool> blnimg(depthMap.Frame()); blnimg.Fill(false);
   for(Array2dIter5C<RealT, bool, bool, bool, bool> it(depthMap, rnimg, bnimg, brnimg, blnimg); it; it++) {

      if(it.Data1()<=0)
         continue;

      // Get neighbours
      Index2dC rn = it.Index()+Index2dC(0,1);
      Index2dC bn = it.Index()+Index2dC(1,0);
      Index2dC brn = it.Index()+Index2dC(1,1);
      Index2dC bln = it.Index()+Index2dC(1,-1);

      // skip if not enough candidates for connected neighbours
      if((!depthMap.Contains(rn) || depthMap[rn]<=0) &&
          (!depthMap.Contains(bn) || depthMap[bn]<=0) &&
          (!depthMap.Contains(brn) || depthMap[brn]<=0) &&
          (!depthMap.Contains(bln) || depthMap[bln]<=0))
         continue;

      // Compute 3D coordinates of current point
      UIntT ind = vlookup[it.Index()];
      Vector3dC pt = v[ind];

      // Process right neighbour
      if(depthMap.Contains(rn) && depthMap[rn]>0) {
         // Compute 3D point coordinates
         ind = vlookup[rn];
         Vector3dC rpt = v[ind];

         // Check connectivity
         it.Data2() = AreConnected(PairC<Vector3dC>(pt, rpt), triMaps, cameras, activeWindow,step);
      }

      // Process bottom neighbour
      if(depthMap.Contains(bn) && depthMap[bn]>0) {
         // Compute 3D point coordinates
         ind = vlookup[bn];
         Vector3dC bpt = v[ind];

         // Check connectivity
         it.Data3() = AreConnected(PairC<Vector3dC>(pt, bpt), triMaps, cameras, activeWindow, step);
      }

      // Process bottom-right neighbour
      if(depthMap.Contains(brn) && depthMap[brn]>0) {
         // Compute 3D point coordinates
         ind = vlookup[brn];
         Vector3dC brpt = v[ind];

         // Check connectivity
         it.Data4() = AreConnected(PairC<Vector3dC>(pt, brpt), triMaps, cameras, activeWindow, step);
      }

      // Process bottom-left neighbour
      if(depthMap.Contains(bln) && depthMap[bln]>0) {
         // Compute 3D point coordinates
         ind = vlookup[bln];
         Vector3dC blpt = v[ind];

         // Check connectivity
         it.Data5() = AreConnected(PairC<Vector3dC>(pt, blpt), triMaps, cameras, activeWindow, step);
      }
   }

   // COMPUTE DEPTH AT PIXEL VERTICES
   ImageC<RealT> tlvimg(depthMap.Frame()); tlvimg.Fill(-1.0);
   ImageC<RealT> trvimg(depthMap.Frame()); trvimg.Fill(-1.0);
   ImageC<RealT> blvimg(depthMap.Frame()); blvimg.Fill(-1.0);
   ImageC<RealT> brvimg(depthMap.Frame()); brvimg.Fill(-1.0);
   for(Array2dIter5C<RealT, RealT, RealT, RealT, RealT> it(depthMap, tlvimg, trvimg, blvimg, brvimg); it; it++) {

      RealT &cpix = it.Data1();
      RealT &tlv = it.Data2();
      RealT &trv = it.Data3();
      RealT &blv = it.Data4();
      RealT &brv = it.Data5();

      if(cpix<=0.0)
         continue;

      tlv = cpix; IntT tlcount = 1;
      trv = cpix; IntT trcount = 1;
      blv = cpix; IntT blcount = 1;
      brv = cpix; IntT brcount = 1;

      Index2dC n;

      // Process top-left pixel neighbour
      n = it.Index()+Index2dC(-1,-1);
      if(brnimg.Contains(n) && brnimg[n]==true) {
         tlv += depthMap[n]; tlcount++;
      }

      // Process top-right pixel neighbour
      n = it.Index()+Index2dC(-1,1);
      if(blnimg.Contains(n) && blnimg[n]==true) {
         trv += depthMap[n]; trcount++;
      }

      // Process bottom-left pixel neighbour
      n = it.Index()+Index2dC(1,-1);
      if(blnimg.Contains(it.Index()) && blnimg[it.Index()]==true) {
         blv += depthMap[n]; blcount++;
      }

      // Process bottom-right pixel neighbour
      n = it.Index()+Index2dC(1,1);
      if(brnimg.Contains(it.Index()) && brnimg[it.Index()]==true) {
         brv += depthMap[n]; brcount++;
      }

      // Process top pixel neighbour
      n = it.Index()+Index2dC(-1,0);
      if(bnimg.Contains(n) && bnimg[n]==true) {
         tlv += depthMap[n]; tlcount++;
         trv += depthMap[n]; trcount++;
      }

      // Process left pixel neighbour
      n = it.Index()+Index2dC(0,-1);
      if(rnimg.Contains(n) && rnimg[n]==true) {
         tlv += depthMap[n]; tlcount++;
         blv += depthMap[n]; blcount++;
      }

      // Process bottom pixel neighbour
      n = it.Index()+Index2dC(1,0);
      if(bnimg.Contains(it.Index()) && bnimg[it.Index()]==true) {
         blv += depthMap[n]; blcount++;
         brv += depthMap[n]; brcount++;
      }

      // Process right pixel neighbour
      n = it.Index()+Index2dC(0,1);
      if(rnimg.Contains(it.Index()) && rnimg[it.Index()]==true) {
         trv += depthMap[n]; trcount++;
         brv += depthMap[n]; brcount++;
      }

      tlv /= tlcount;
      trv /= trcount;
      blv /= blcount;
      brv /= brcount;
   }

   UIntT nv2;
   DListC<Vector3dC> v2;
   nv2 = 0;

   // Process each top-left pixel vertex
   ImageC<UIntT> tlvlookup(depthMap.Frame());
   tlvlookup.Fill(0);
   for(Array2dIter2C<RealT, UIntT> it(tlvimg, tlvlookup); it; it++) {

      if(it.Data1()>0) {
         // Get image point coordinates
         Vector2dC imgPt(it.Index().Col()-0.5, it.Index().Row()-0.5);

         // Determine the 3D direction of the ray
         Vector3dC dir;
         refCam.ProjectInverseDirection(dir, imgPt);
         dir.MakeUnit();

         // Compute 3D point coordinates
         Vector3dC pt = camOrg + dir*it.Data1();

         // Store data
         v2.InsLast(pt);
         it.Data2() = nv2;
         nv2++;
      }
   }

   // Process each top-right pixel vertex
   ImageC<UIntT> trvlookup(depthMap.Frame());
   trvlookup.Fill(0);
   for(Array2dIter2C<RealT, UIntT> it(trvimg, trvlookup); it; it++) {

      if(it.Data1()>0) {
         Index2dC n = it.Index()+Index2dC(0,1);
         if(tlvimg.Contains(n) && tlvimg[n]==it.Data1()) {
            // The point already exists

            // Access previous point and store data
            it.Data2() = tlvlookup[n];

            continue;
         }

         // A new point must be created

         // Get image point coordinates
         Vector2dC imgPt(it.Index().Col()+0.5, it.Index().Row()-0.5);

         // Determine the 3D direction of the ray
         Vector3dC dir;
         refCam.ProjectInverseDirection(dir, imgPt);
         dir.MakeUnit();

         // Compute 3D point coordinates
         Vector3dC pt = camOrg + dir*it.Data1();

         // Store data
         v2.InsLast(pt);
         it.Data2() = nv2;
         nv2++;
      }
   }

   // Process each bottom-left pixel vertex
   ImageC<UIntT> blvlookup(depthMap.Frame());
   blvlookup.Fill(0);
   for(Array2dIter2C<RealT, UIntT> it(blvimg, blvlookup); it; it++) {

      if(it.Data1()>0) {
         Index2dC n = it.Index()+Index2dC(1,0);
         if(tlvimg.Contains(n) && tlvimg[n]==it.Data1()) {
            // The point already exists

            // Access previous point and store data
            it.Data2() = tlvlookup[n];

            continue;
         }

         n = it.Index()+Index2dC(1,-1);
         if(trvimg.Contains(n) && trvimg[n]==it.Data1()) {
            // The point already exists

            // Access previous point and store data
            it.Data2() = trvlookup[n];

            continue;
         }

         // A new point must be created

         // Get image point coordinates
         Vector2dC imgPt(it.Index().Col()-0.5, it.Index().Row()+0.5);

         // Determine the 3D direction of the ray
         Vector3dC dir;
         refCam.ProjectInverseDirection(dir, imgPt);
         dir.MakeUnit();

         // Compute 3D point coordinates
         Vector3dC pt = camOrg + dir*it.Data1();

         // Store data
         v2.InsLast(pt);
         it.Data2() = nv2;
         nv2++;
      }
   }

   // Process each bottom-right pixel vertex
   ImageC<UIntT> brvlookup(depthMap.Frame());
   brvlookup.Fill(0);
   for(Array2dIter2C<RealT, UIntT> it(brvimg, brvlookup); it; it++) {

      if(it.Data1()>0) {
         Index2dC n = it.Index()+Index2dC(1,1);
         if(tlvimg.Contains(n) && tlvimg[n]==it.Data1()) {
            // The point already exists

            // Access previous point and store data
            it.Data2() = tlvlookup[n];

            continue;
         }

         n = it.Index()+Index2dC(1,0);
         if(trvimg.Contains(n) && trvimg[n]==it.Data1()) {
            // The point already exists

            // Access previous point and store data
            it.Data2() = trvlookup[n];

            continue;
         }

         n = it.Index()+Index2dC(0,1);
         if(blvimg.Contains(n) && blvimg[n]==it.Data1()) {
            // The point already exists

            // Access previous point and store data
            it.Data2() = blvlookup[n];

            continue;
         }

         // A new point must be created

         // Get image point coordinates
         Vector2dC imgPt(it.Index().Col()+0.5, it.Index().Row()+0.5);

         // Determine the 3D direction of the ray
         Vector3dC dir;
         refCam.ProjectInverseDirection(dir, imgPt);
         dir.MakeUnit();

         // Compute 3D point coordinates
         Vector3dC pt = camOrg + dir*it.Data1();

         // Store data
         v2.InsLast(pt);
         it.Data2() = nv2;
         nv2++;
      }
   }

   DListC<UIntT> faceInd;
   DListC<ByteRGBAValueC> colours;
   for(Array2dIterC<RealT> it(depthMap); it; it++) {

      if(it.Data()>0) {
         UIntT ntl= tlvlookup[it.Index()];
         UIntT ntr= trvlookup[it.Index()];
         UIntT nbl= blvlookup[it.Index()];
         UIntT nbr= brvlookup[it.Index()];

         // insert 2 triangles
         faceInd.InsLast(ntr); faceInd.InsLast(ntl); faceInd.InsLast(nbr);
         faceInd.InsLast(ntl); faceInd.InsLast(nbl); faceInd.InsLast(nbr);

         if(embedColour==true) {
            // insert colour (same colour for both triangles) and alpha
            ByteRGBAValueC ncol(fgRefImage[it.Index()].Red(), fgRefImage[it.Index()].Green(),
                                fgRefImage[it.Index()].Blue(), alphaMatte.IsEmpty()?255:alphaMatte[it.Index()]);
            colours.InsLast(ncol);
            colours.InsLast(ncol);
         }
      }
   }

   // Save mesh
   if(!output.IsEmpty()) {
      char pFile[255];
      sprintf(pFile,output.chars(),refCamID,frameNumber);
      OStreamC meshFile(pFile);

      if(embedColour==true) {
         meshFile << 'C' << '\n';
      }
      meshFile << v2.Size() << " " << faceInd.Size()/3 << '\n';
      for(DLIterC<Vector3dC> it(v2); it; it++) {
         meshFile << it.Data() << '\n';
      }
      for(DLIterC<UIntT> it(faceInd); it; it++) {
         if(embedColour==true) {
            meshFile << 7;
         } else {
            meshFile << 3;
         }
         meshFile << " " << it.Data();
         it++;
         meshFile << " " << it.Data();
         it++;
         meshFile << " " << it.Data();
         if(embedColour==true) {
            ByteRGBAValueC colour = colours.PopFirst();
            meshFile << " " << (int) colour.Red() << " " << (int) colour.Green() << " "
                     << (int) colour.Blue() << " " << (int) colour.Alpha();
         }
         meshFile << '\n';
      }
      meshFile.Close();
   }

   return true;
}

bool
AreConnected(const PairC<Vector3dC>& pts,
             const SArray1dC<ImageC<ByteT> >& triMaps,
             const PinholeCameraArrayC& cameras,
             const IndexRange2dC activeWindow,
             RealT step)
{

   // Compute distance between the two points
   RealT dist_max = pts.Data1().EuclidDistance(pts.Data2());

   bool res = true;

   // Compute unit vector on line connecting the two points
   Vector3dC dir = pts.Data2()-pts.Data1();
   dir.MakeUnit();

   for(RealT dist=0+step; dist<=dist_max-step; dist+=step) {
      Vector3dC pt = pts.Data1() + dir*dist;
      // Check if point stays within the visual hull
      if(!IsWithinVisualHull(pt, cameras, triMaps, activeWindow, false)) {
         res = false;
         break;
      }
   }

   /*if(dist_max>5*step) {
      res = false;
   }*/

   return res;
}

// Interpolate colour from four pixel neighbours
bool
InterpolateColor(const ImageC<ByteRGBValueC>& rgbImage,
                 const Vector2dC& imgPt,
                 RealRGBValueC& colour)
{
   // Bi-linear interpolation
   IntT fr = Floor(imgPt.Y());
   IntT fc = Floor(imgPt.X());
   RealT u = imgPt.Y() - fr;
   RealT t = imgPt.X() - fc;
   RealT onemu = (1.0-u);
   RealT onemt = (1.0-t);
   const ByteRGBValueC* pixel1 = &rgbImage[fr][fc];
   const ByteRGBValueC* pixel2 = &rgbImage[fr+1][fc];

   ByteRGBValueC p10 = pixel1[0];
   ByteRGBValueC p11 = pixel1[1];
   ByteRGBValueC p20 = pixel2[0];
   ByteRGBValueC p21 = pixel2[1];

   colour = (p10*onemt*onemu)+(p11*t*onemu)+(p20*onemt*u)+(p21*t*u);

   return true;
}

// Project auxiliary cameras to the reference camera
bool
ProjectToReference(const SArray1dC<ImageC<ByteRGBValueC> >& rgbImages,
                   const ImageC<RealT>& depthMap,
                   const PinholeCameraArrayC& cameras,
                   IntT refCamID,
                   SArray1dC<IntT>& auxCamIDs,
                   SArray1dC<ImageC<RealRGBValueC> >& warpedImages,
                   SArray1dC<ImageC<RealT> >& errorMaps,
                   const IndexRange2dC activeWindow,
                   IntT frameNumber,
                   StringC outputSTATS,
                   bool useRobustStats,
                   bool useIsotropic,
                   bool VERBOSE)
{

   warpedImages = SArray1dC<ImageC<RealRGBValueC> >(rgbImages.Size());
   errorMaps = SArray1dC<ImageC<RealT> >(rgbImages.Size());
   for(Array1dIterC<IntT> it(auxCamIDs); it; it++) {
      warpedImages[it.Data()] = ImageC<RealRGBValueC>(depthMap.Frame());
      warpedImages[it.Data()].Fill(RealRGBValueC(0.0, 0.0, 0.0));
      errorMaps[it.Data()] = ImageC<RealT>(depthMap.Frame());
      errorMaps[it.Data()].Fill(0.0);
   }
   // Get reference camera
   const PinholeCameraC& refCam = cameras[refCamID];
   Vector3dC refCamOrg;
   refCam.Origin(refCamOrg);

   SArray1dC<DListC<VectorC> > samples(auxCamIDs.Size());

   char pFile[255];

   // Iterate over the depth map
   for (Array2dIterC<RealT> itimg(depthMap); itimg; itimg++)
   {
      if(itimg.Data()>0.0) {
         // Get depth of image point
         RealT& depth = itimg.Data();
         Vector2dC imgPt(itimg.Index().Col(), itimg.Index().Row());

         // Determine the 3D direction of the ray
         Vector3dC dir;
         refCam.ProjectInverseDirection(dir, imgPt);
         dir.MakeUnit();

         // Compute coordinates of 3D point
         Vector3dC Pt = refCamOrg + dir*depth;

         // project into auxiliary camera
         SArray1dC<Vector2dC> imgPts;
         ComputeProjections(imgPts, Pt, cameras, auxCamIDs);

         for(Array1dIter3C<IntT, Vector2dC, DListC<VectorC> > it(auxCamIDs, imgPts, samples); it; it++) {
            if (it.Data2().X()>activeWindow.LCol() && it.Data2().X()<activeWindow.RCol() && it.Data2().Y()>activeWindow.TRow() && it.Data2().Y()<activeWindow.BRow()) {
               RealRGBValueC auxC;
               InterpolateColor(rgbImages[it.Data1()], it.Data2(), auxC);
               warpedImages[it.Data1()][itimg.Index()] = auxC;
               ByteRGBValueC refC = rgbImages[refCamID][itimg.Index()];
               errorMaps[it.Data1()][itimg.Index()] = Sqrt((Sqr(refC.Red()-auxC.Red())+Sqr(refC.Green()-auxC.Green())+Sqr(refC.Blue()-auxC.Blue()))/3.0);
               it.Data3().InsLast(VectorC(refC.Red()-auxC.Red(), refC.Green()-auxC.Green(), refC.Blue()-auxC.Blue()));
           }
         }

      }
   }

   if(VERBOSE) {
      for(Array1dIterC<IntT> it(auxCamIDs); it; it++) {
         Save("@X:warpedImage"+StringC(it.Data()),warpedImages[it.Data()]);
         Save("@X:errorMap"+StringC(it.Data()),errorMaps[it.Data()]);
      }
   }

   for(Array1dIter2C<IntT, DListC<VectorC> > it(auxCamIDs, samples); it; it++) {
      VectorC mean;
      MatrixRSC cov;
      if(useRobustStats) {
         // compute mean
         RobustMeanCovarianceC(it.Data2(), mean, cov, true);
      } else {
         // compute mean
         MeanCovarianceC pixelStats(it.Data2(), true);
         mean = pixelStats.Mean();
         DListC<VectorC> zlist;
         // force mean to zero
         for(DLIterC<VectorC> itl(it.Data2()); itl; itl++) {
            zlist.InsLast(itl.Data()-mean);
         }
         // re-compute statistics with zero mean
         pixelStats = MeanCovarianceC(zlist, true);
         mean = pixelStats.Mean();
         cov = pixelStats.Covariance();
      }
      if(useIsotropic) {
         RealT var = (cov[0][0]+cov[1][1]+cov[2][2])/3.0;
         cov[0][0] = var;
         cov[1][1] = var;
         cov[2][2] = var;
      }

      Vector3dC meanStd(Sqrt(cov[0][0]), Sqrt(cov[1][1]), Sqrt(cov[2][2]));
      sprintf(pFile,outputSTATS.chars(),refCamID,it.Data1(),frameNumber);
      OStreamC fout(pFile);
      fout << meanStd << endl;
      fout.Close();

      if(VERBOSE) {
         cerr << "Statistics for camera " << it.Data1() << " to " << refCamID << endl;
         cerr << "Per channel standard deviation: " << meanStd << endl;
      }
   }

   return true;

}

bool RobustMeanCovarianceC(const DListC<VectorC> & data,
                           VectorC& mean,
                           MatrixRSC& cov,
                           bool zeroMean) 
{
   IntT number = data.Size();
   if(number == 0)
      return false;

   if(!zeroMean) {
      CollectionC<OrderedVectorC> smean(number);
      for(DLIterC<VectorC> it(data); it; it++) {
         smean.Insert(it.Data());
      }
      mean = smean.KthHighest(number/2);
   } else {
      mean = VectorC(0.0,0.0,0.0);
   }

   SArray2dC<CollectionC<RealT> > scov(3,3);
   scov[0][0] = CollectionC<RealT>(number);
   scov[0][1] = CollectionC<RealT>(number);
   scov[0][2] = CollectionC<RealT>(number);
   scov[1][1] = CollectionC<RealT>(number);
   scov[1][2] = CollectionC<RealT>(number);
   scov[2][2] = CollectionC<RealT>(number);

   for(DLIterC<VectorC> it(data); it; it++) {
      VectorC diff = mean-it.Data();
      MatrixRSC diffmat= diff.OuterProduct();

      scov[0][0].Insert(diffmat[0][0]);
      scov[0][1].Insert(diffmat[0][1]);
      scov[0][2].Insert(diffmat[0][2]);
      scov[1][1].Insert(diffmat[1][1]);
      scov[1][2].Insert(diffmat[1][2]);
      scov[2][2].Insert(diffmat[2][2]);
   }

   SArray2dC<RealT> mcov(3,3);
   mcov[0][0] = scov[0][0].KthHighest(number/2);
   mcov[0][1] = scov[0][1].KthHighest(number/2);
   mcov[0][2] = scov[0][2].KthHighest(number/2);
   mcov[1][1] = scov[1][1].KthHighest(number/2);
   mcov[1][2] = scov[1][2].KthHighest(number/2);
   mcov[2][2] = scov[2][2].KthHighest(number/2);
 
   cov = MatrixC(mcov[0][0], mcov[0][1], mcov[0][2],
                 mcov[0][1], mcov[1][1], mcov[1][2],
                 mcov[0][2], mcov[1][2], mcov[2][2]);

   return true;
}


