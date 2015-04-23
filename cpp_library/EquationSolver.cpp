/**
 * File:   EquationSolver.cpp
 * Author: kenobi
 *
 * Created on April 22, 2015, 12:51 PM
 */

#include "EquationSolver.hpp"
#include <cmath>

namespace nde {

EquationSolver::EquationSolver() {
	;
}


double EquationSolver::secant(double x0, double x1) {

	double y0 = this->solverF(x0);
	double y1 = this->solverF(x1);

	// the solution is stored in the (x,y) values
	double y;
	double x;

	if (std::abs(y0) < std::abs(y1)) {
		y = y0;
		x = x0;
	} else {
		y = y1;
		x = x1;
	}

	while (std::abs(y) > 1e-12) {
		x = x0 - (x1 - x0) / (y1 / y0 - 1.0);
		y = this->solverF(x);
		if (std::abs(y0) < std::abs(y1)) {
			x1 = x;
			y1 = y;
		} else {
			x0 = x;
			y0 = y;
		}
	}

	return x;

}

} /* end of namespace nde */

