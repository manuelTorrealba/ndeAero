/**
 * File:   TestODESolver.cpp
 * Author: kenobi
 *
 * Created on May 7, 2015, 9:53 AM
 */

#include "TestODESolver.hpp"

namespace nde {

SimpleExponential::SimpleExponential(double lambda, unsigned int order) :
						ODESolver(order, 1.), _lambda(lambda) {
	;
}

Vector<double> SimpleExponential::odeSolverDy(double t,
														const Vector<double>& y) const {
	return Vector<double>(y.size(), -y(0) * _lambda);
}

void testODESolver() {

	// initial solution
	Vector<double> y0(1, 1.0);

	std::cout << "EULER" << std::endl;
	SimpleExponential simple_exponential_euler(1.0, 1);
	Matrix<double> y_euler = simple_exponential_euler.solve(0., 1., 0.1, y0);
	for (unsigned int i = 0; i < y_euler.numcols(); ++i)
		std::cout << y_euler(0,i) << std::endl;

	std::cout << "RK4" << std::endl;
	SimpleExponential simple_exponential_rk4(1.0, 4);
	Matrix<double> y_rk4 = simple_exponential_rk4.solve(0., 1., 0.1, y0);
	for (unsigned int i = 0; i < y_rk4.numcols(); ++i)
		std::cout << y_rk4(0,i) << std::endl;

	std::cout << "RK5" << std::endl;
	SimpleExponential simple_exponential_rk5(1.0, 6);
	Matrix<double> y_rk5 = simple_exponential_rk5.solve(0., 1., 0.1, y0);
	for (unsigned int i = 0; i < y_rk5.numcols(); ++i)
		std::cout << y_rk5(0,i) << std::endl;

	return;

}

} /* end of namespace nde */

