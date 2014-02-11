# defs.mk Lip Tracking File
USESLIBS = RavlCore RavlIO RavlImageIO RavlMath RavlVideoIO RavlDPDisplay RavlImageProc RavlExtImgIO RobustClustering Bspline RavlPatternRec 

PROGLIBS = OmniLib SegLib ParticleFilter vbapAT
HEADERS = EigenLoadParams.hh
#MAINS = FOPropJMeasTracker.cc OFPropLLmeasTracker.cc AREstimationUsingLandmarks.cc SOARPropLLMeasTracker.cc
#MAINS = SOARPropLLMeasTracker.cc
#MAINS = EJTestSetup.cc
#MAINS = FOPropJMeasTracker.cc
#MAINS = BootstrapProcess.cc
MAINS = FOJSRTracker.cc OFLLSRTracker.cc FOEigenLLTracker.cc
#OFLLSRTracker.cc StateVectorTester.cc
