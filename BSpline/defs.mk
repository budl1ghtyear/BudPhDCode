 # defs.mk
 PLIB = Bspline 
 USESLIBS = RavlCore RavlIO RavlImageIO RavlVideoIO RavlDPDisplay RavlImageProc RavlExtImgIO 
 HEADERS = BSplineC.hh BasisSpline.hh MBspline.hh
 #BasisSpline.hh URSBSpline.hh 
 SOURCES = BSplineC.cc BasisSpline.cc MBspline.cc
 #BasisSpline.cc URSBSpline.cc
 #MAINS = FitClosedSpline.cc MBSplineFitter.cc
 #MAINS = MBSplineFitter.cc
 MAINS = TestSplineClasses.cc
