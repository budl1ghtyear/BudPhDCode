# defs.mk
#PARTICLE FILTERING FRAMEWORK
#PLIB = ParticleFilter
#USESLIBS = RavlCore RavlIO RavlImageIO RavlVideoIO RavlImage Bspline Modelling RavlOpenCV RavlDPDisplay
USESLIBS= RavlVideoIO RavlDPDisplay RavlImageIO RavlExtImgIO RavlImage RavlCore RavlImageIO RavlImageProc Bspline Modelling RavlLibFFmpeg RavlOpenCV
#RavlOpenCV 
HEADERS = BaseStateVector.hh Particle.hh BasePropagationModel.hh BaseMeasurementModel.hh BaseResamplingModel.hh BaseParticleFilter.hh SIR.hh SystematicResampling.hh PCAStateVector.hh ZOAVPropModel.hh FOPropagationModel.hh ICAStateVector.hh EdgeMeasurementModel.hh BANCAGT.hh TrackingUtilityFunctions.hh GroupSpline.hh HistMeasurementModel.hh ChiHistMeasurementModel.hh BoxStateVector.hh BoxPropModel.hh
SOURCES = BaseStateVector.cc Particle.cc SIR.cc SystematicResampling.cc PCAStateVector.cc ZOAVPropModel.cc FOPropagationModel.cc ICAStateVector.cc EdgeMeasurementModel.cc BANCAGT.cc TrackingUtilityFunctions.cc GroupSpline.cc HistMeasurementModel.cc ChiHistMeasurementModel.cc BoxStateVector.cc BoxPropModel.cc
MAINS=GroupSplineTest.cc
#MAINS=loadvideo.cc myconv.cc
#HEADERS = BaseStateVector.hh Particle.hh BasePropagationModel.hh BaseMeasurementModel.hh BaseResamplingModel.hh BaseParticleFilter.hh SIR.hh PointStateVector.hh EigenStateVector.hh SOARPropagationModel.hh LLMeasurementModel.hh SystematicResampling.hh TrackingFunctions.hh
#SOURCES = BaseStateVector.cc Particle.cc SIR.cc PointStateVector.cc EigenStateVector.cc SOARPropagationModel.cc LLMeasurementModel.cc SystematicResampling.cc TrackingFunctions.cc 
#HEADERS = BaseStateVector.hh Particle.hh EigenStateVector.hh
#SOURCES = BaseStateVector.cc Particle.cc EigenStateVector.cc
