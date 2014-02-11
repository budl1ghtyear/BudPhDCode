PACKAGE = Ravl/Image

REQUIRES = OmniLib

HEADERS = FaceLocate.hh

USESLIBS = OmniLib RavlSVM RavlImageProc RavlImage RavlCore

ifdef SHAREDBUILD

  SOURCES = FaceLocate.cc

  #TESTEXES = testFaceLocate.cc

  #EXAMPLES = exFaceLocate.cc
	
  MAINS = testFaceLocate.cc exFaceLocate.cc

  PLIB = RavlOmniProc

  PROGLIBS = RavlExtImgIO

  AUXFILES = face.jpg

  AUXDIR = share/RAVL/testData

  CCPPFLAGS += -DOMNIHOME=\"$(OMNIHOME)\"

endif
