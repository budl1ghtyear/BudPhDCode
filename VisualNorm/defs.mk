PACKAGE =Vtefm 
USESLIBS=RavlCore RavlIO Ravl3D RavlImage RavlDPDisplay Ravl3DIO RavlImageIO RavlExtImgIO RavlGUI3D RavlMath RavlPatternRec RavlImageProc Auto SegLib
MAINS = testColorTransformGuided.cc testStirlingData.cc testHistogramSpecification.cc
SOURCES = ColorNorm.cc ColorNormHist.cc ColorSpaces.cc
HEADERS = ColorNorm.hh ColorSpaces.hh
PLIB=VtefmVisualNorm


