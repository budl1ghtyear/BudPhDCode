PLIB = SegLib
SOURCES = MCD.cc CascadedMCD.cc KSmirnovMCD.cc SpectralClustering.cc CostFunctions.cc MouthRegion.cc ColourConvert.cc SampleData.cc XMLToGTConv.cc LipSegmentation.cc BANCAGT.cc LipExtractorC.cc DesignFuzzyCMeansCluster.cc
HEADERS = MCD.hh CascadedMCD.hh KSmirnovMCD.hh SpectralClustering.hh CostFunctions.hh MouthRegion.hh LipClusteringC.hh ColourConvert.hh SampleData.hh XMLToGTConv.hh LipSegmentation.hh ColourSpaceImage.hh BANCAGT.hh LipExtractorC.hh DesignFuzzyCMeansCluster.hh
USESLIBS = RavlCore RavlIO RavlImage RavlDPDisplay RavlImageIO  RavlExtImgIO RavlMath RavlMath RavlPatternRec RavlImageProc Auto ChiSquared RavlOpenCV
#PROGLIBS = RavlOpenCV
#CCPPFLAGS += -DOMNIHOME=\"$(OMNIHOME)\"
