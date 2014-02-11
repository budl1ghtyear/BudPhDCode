# defs.mk
#PARTICLE FILTERING FRAMEWORK
PLIB = ParticleFilter
USESLIBS = RavlCore RavlIO RavlImageIO RavlVideoIO RavlImage Bspline Modelling RavlOpenCV SegLib
HEADERS = BaseStateVector.hh Particle.hh BasePropagationModel.hh BaseMeasurementModel.hh BaseResamplingModel.hh BaseParticleFilter.hh SIR.hh SystematicResampling.hh PCAStateVector.hh ZOAVPropModel.hh FOPropagationModel.hh ICAStateVector.hh ACMeasurementModel.hh TrackingUtilityFunctions.hh EdgeMeasurementModel.hh
SOURCES = BaseStateVector.cc Particle.cc SIR.cc SystematicResampling.cc PCAStateVector.cc ZOAVPropModel.cc FOPropagationModel.cc ICAStateVector.cc ACMeasurementModel.cc TrackingUtilityFunctions.cc EdgeMeasurementModel.cc
MAINS = TestPCAStateVector.cc ZOAVACPCATracker.cc 
#HEADERS = BaseStateVector.hh Particle.hh BasePropagationModel.hh BaseMeasurementModel.hh BaseResamplingModel.hh BaseParticleFilter.hh SIR.hh PointStateVector.hh EigenStateVector.hh SOARPropagationModel.hh LLMeasurementModel.hh SystematicResampling.hh TrackingFunctions.hh
#SOURCES = BaseStateVector.cc Particle.cc SIR.cc PointStateVector.cc EigenStateVector.cc SOARPropagationModel.cc LLMeasurementModel.cc SystematicResampling.cc TrackingFunctions.cc 
#HEADERS = BaseStateVector.hh Particle.hh EigenStateVector.hh
#SOURCES = BaseStateVector.cc Particle.cc EigenStateVector.cc
