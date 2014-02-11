# defs.mk
#PARTICLE FILTERING FRAMEWORK
PLIB = ParticleFilter
USESLIBS = RavlCore RavlIO RavlImageIO RavlVideoIO RavlImage Bspline Modelling
HEADERS = BaseStateVector.hh Particle.hh BasePropagationModel.hh BaseMeasurementModel.hh BaseResamplingModel.hh BaseParticleFilter.hh SIR.hh SystematicResampling.hh PCAStateVector.hh ZOAVPropModel.hh FOPropagationModel.hh
SOURCES = BaseStateVector.cc Particle.cc SIR.cc SystematicResampling.cc PCAStateVector.cc ZOAVPropModel.cc FOPropagationModel.cc
MAINS = TestPCAStateVector.cc

