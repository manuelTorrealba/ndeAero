/**
 * File:   TestThinAirfoil.cpp
 * Author: kenobi
 *
 * Created on April 23, 2015, 09:10 AM
 */

#include "TestThinAirfoil.hpp"

namespace nde {

ParabolicArc::ParabolicArc(unsigned int n, double curvature)
								:ThinAirfoil(n), _curvature(curvature) {
	;
}

double ParabolicArc::camber(double t) const {
	return 4.0 * _curvature * (1.0 - t) * t;
}

double ParabolicArc::dCamberDx(double t) const {
	return 4.0 * _curvature * (1.0 - 2.0 * t);
}

} /* end of namespace nde */

