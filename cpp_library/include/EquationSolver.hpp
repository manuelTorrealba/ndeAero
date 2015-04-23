/**
 * File:   EquationSolver.hpp
 * Author: kenobi
 *
 * Created on April 22, 2015, 12:51 PM
 */

#ifndef INCLUDE_EQUATION_SOLVER_HPP
#define INCLUDE_EQUATION_SOLVER_HPP

namespace nde {

class EquationSolver {
public:
	EquationSolver();
	virtual	double solverF(double x) = 0;

protected:
	double secant(double x0, double x1);

};

} /* end of namespace nde */

#endif

