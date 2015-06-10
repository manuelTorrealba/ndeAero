/**
 * File:   TestODESolver.cpp
 * Author: kenobi
 *
 * Created on May 7, 2015, 9:53 AM
 */

#include "Test.hpp"
#include "TestODESolver.hpp"

#include <iostream>

namespace nde {

/**
  * SimpleExponential class
  */
SimpleExponential::SimpleExponential(double lambda, unsigned int order) :
						ODESolver(order, 1.), _lambda(lambda) {
	;
}

Vector<double> SimpleExponential::odeSolverDy(double t,
														const Vector<double>& y) const {
	return Vector<double>(y.size(), -y(0) * _lambda);
}


/******************************************************************************/
/* Blasius boundary Layer
/******************************************************************************/

Blasius::Blasius(unsigned int order) : ODESolver(order, 1.) {
	;
}

Vector<double> Blasius::odeSolverDy(double t, const Vector<double>& y) const {

	Vector<double> dy(3);
	dy(0) = y(1);
	dy(1) = y(2);
	dy(2) = - 0.5 * y(0) * y(2);

	return dy;

};




bool testODESolver() {

	std::cout << "/******************************************************/" << std::endl;
	std::cout << "/* ODE Solver test" << std::endl;
	std::cout << "/******************************************************/" << std::endl;

	{

	std::cout << "SIMPLE EXPONENTIAL dY = - Y dt" << std::endl;

	// initial solution
	Vector<double> y0(1, 1.0);

	std::cout << "EULER" << std::endl;
	SimpleExponential simple_exponential_euler(1.0, 1);
	Matrix<double> y_euler = simple_exponential_euler.solve(0., 1., 10, y0);
	for (unsigned int i = 0; i < y_euler.numcols(); ++i)
		std::cout << y_euler(0,i) << std::endl;

	std::cout << "RK4" << std::endl;
	SimpleExponential simple_exponential_rk4(1.0, 4);
	Matrix<double> y_rk4 = simple_exponential_rk4.solve(0., 1., 10, y0);
	for (unsigned int i = 0; i < y_rk4.numcols(); ++i)
		std::cout << y_rk4(0,i) << std::endl;

	std::cout << "RK5" << std::endl;
	SimpleExponential simple_exponential_rk5(1.0, 6);
	Matrix<double> y_rk5 = simple_exponential_rk5.solve(0., 1., 10, y0);
	for (unsigned int i = 0; i < y_rk5.numcols(); ++i)
		std::cout << y_rk5(0,i) << std::endl;

	}


	{

	std::cout << "BLASIUS BOUNDARY LATER 1/2 f f'' + f''' = 0" << std::endl;

	// initial solution
	Vector<double> y0(3);
	y0(0) = 0.;
	y0(1) = 0.;
	y0(2) = 0.3321; // taken from some text

	std::cout << "EULER" << std::endl;
	Blasius blasius_euler(1);
	Matrix<double> y_euler = blasius_euler.solve(0., 8., 20, y0);
	for (unsigned int i = 0; i < y_euler.numcols(); ++i)
		std::cout << y_euler(0,i) << "," << y_euler(1,i)
					 << "," << y_euler(2,i) << std::endl;

	std::cout << "RK4" << std::endl;
	Blasius blasius_rk4(4);
	Matrix<double> y_rk4 = blasius_rk4.solve(0., 8., 20, y0);
	for (unsigned int i = 0; i < y_rk4.numcols(); ++i)
		std::cout << y_rk4(0,i) << "," << y_rk4(1,i)
					 << "," << y_rk4(2,i) << std::endl;

	std::cout << "RK5" << std::endl;
	Blasius blasius_rk5(6);
	Matrix<double> y_rk5 = blasius_rk5.solve(0., 8., 20, y0);
	for (unsigned int i = 0; i < y_rk5.numcols(); ++i)
		std::cout << y_rk5(0,i) << "," << y_rk5(1,i)
					 << "," << y_rk5(2,i) << std::endl;

	}

	return true;

}

} /* end of namespace nde */

