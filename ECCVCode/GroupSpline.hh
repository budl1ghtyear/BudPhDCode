#ifndef GROUPSPLINE_HH
#define GROUPSPLINE_HH

//Define a set of methods required to do group spline based particle filter stuff
#include "Ravl/Array1d.hh"
#include "Ravl/SArray1d.hh"
#include "Particle.hh"
#include "Ravl/PatternRec/Sample.hh"
#include "Ravl/PatternRec/SampleIter.hh"
#include "Ravl/LinePP2d.hh"
#include "Ravl/String.hh"
#include "Ravl/StringList.hh"
#include "Ravl/Vector2d.hh"
#include "Ravl/Point2d.hh"
#include <string>
#include <iostream>
#include <fstream>
namespace RavlN
{
void GroupBSpline(Array1dC<ParticleC> &myparticles, SampleC<SArray1dC<Point2dC> > &curvepts, SampleC<SArray1dC<LinePP2dC> > &normals);
}


#endif
