/**
 * File:   TestThinAirfoilTheory.cpp
 * Author: kenobi
 *
 * Created on Jun 8, 2015, 11:00 PM
 */

#include "Constants.hpp"
#include "TestThinAirfoil.hpp"
#include "ThinAirfoil.hpp"
#include <cmath>
#include <iostream>

namespace nde {

bool testThinAirfoilTheory() {

	std::cout << "/******************************************************/" << std::endl;
	std::cout << "/* Thin Airfoil theory test" << std::endl;
	std::cout << "/******************************************************/" << std::endl;

	nde::ParabolicArc pa(14, 0.1); // 2^14 = 16384 integration points
											// and 0.1 for maximum camber

	double cl_num = pa.calcCL(0.0);
	double cm_num = pa.calcCM(0.0, 0.25);
	std::cout << "C_L = " << cl_num << std::endl;
	std::cout << "C_M 1/4 = " << cm_num << std::endl;

	double cl_analytic = 4.0 * M_PI * 0.1;
	double cm_analytic = -M_PI * 0.1;

	std::cout << "C_L (theory) = " << cl_analytic << std::endl;
	std::cout << "C_M 1/4 (theory) = " << cm_analytic << std::endl;

	if ( (std::abs(cl_num - cl_analytic) < TEST_BIG_TOLERANCE) &&
		  (std::abs(cm_num - cm_analytic) < TEST_BIG_TOLERANCE) )
		return true;
	else
		return false;

}

} /* end of namespace nde */

