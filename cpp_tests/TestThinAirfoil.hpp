/**
 * File:   TestThinAirfoil.hpp
 * Author: kenobi
 *
 * Created on April 23, 2015, 09:10 AM
 */

#ifndef INCLUDE_TEST_THIN_AIRFOIL_HPP
#define INCLUDE_TEST_THIN_AIRFOIL_HPP

#include "ThinAirfoil.hpp"

namespace nde {

class ParabolicArc: public ThinAirfoil {
public:
	ParabolicArc(unsigned int n, double curvature);
	virtual double camber(double t) const;
	virtual	double dCamberDx(double t) const;

private:
	double _curvature;

};

} /* end of namespace nde */

#endif

