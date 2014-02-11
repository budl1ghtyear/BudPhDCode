# defs.mk
#FACTORED SAMPLING 
PLIB = ParticleFilter
USESLIBS = RavlCore RavlIO RavlImageIO RavlVideoIO RavlImage 
HEADERS = BaseStateVector.hh Particle.hh BasePropagationModel.hh BaseMeasurementModel.hh BaseResamplingModel.hh BaseParticleFilter.hh SIR.hh PointStateVector.hh JMeasurementModel.hh FOPropagationModel.hh SystematicResampling.hh
SOURCES = Particle.cc SIR.cc PointStateVector.cc JMeasurementModel.cc FOPropagationModel.cc SystematicResampling.cc

