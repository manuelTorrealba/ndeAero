/**
 * File:   ODESolver.hpp
 * Author: kenobi
 *
 * Created on May 4, 2015, 11:00 PM
 */

#ifndef INCLUDE_ODESOLVER_HPP
#define INCLUDE_ODESOLVER_HPP

#include "Vector.hpp"
#include "Matrix.hpp"

namespace nde {


/**
  * 1D ODE solver
  * dy/dx = f(x,y)
  *
  * Current available methods:
  *
  * Euler: y(n+1) = y(n) + h * f(x_n, y_n)
  *
  * Runge-Kutta (N = 4, 5 and 6 order):
  *
  * y(n+1) = y(n) + sum_(i<=N)(b_i * k_i)
  * k_i = f(x_n + c_i * h, sum_(j<i)(a_ij * k_j))
  *
  */

class ODESolver {
public:
   /**
     * Constructor:
     *
     * n = 1, Euler; 4, RK4; 6, RK6.
     * max_error: maximum error allowed for time steppinf from x to x+h.
     *
     */
	ODESolver(unsigned int n, double max_error);
	
   /**
     * Output:
     * y(0,i): first row is the x-grid.
	  * y(1 + j, i): solution variables in the same order than y0.
     *
     */
	Matrix<double> solve(double x0, double xn, unsigned int n_steps,
							const Vector<double>& y0) const;

	/**
     * f(x,y) function
     */
	virtual Vector<double> odeSolverDy(double x,
												const Vector<double>& y) const = 0;

   /**
     * possibility of introducing a discrete jump at point x.
	  * y+(x+) = y-(x-) + Jump(x-,y-)
     */
	virtual Vector<double> odeSolverJumpy(double x,
													const Vector<double>& y) const = 0;

private:
	unsigned int _n;
	double _max_error;

	Vector<double> nextStep(double x, const Vector<double>& y, double h,
								double& err_estimate) const;

	Vector<double> RungeKutta(unsigned int n, unsigned int order, double x,
									const Vector<double>& y, double h) const;

	/* routine with hard-coded Runge Kutta coefficients	*/
	void RungeKuttaCoefficients(unsigned int n, unsigned int order,
		Vector<double>& a, Matrix<double>& b, Vector<double>& c) const;

};


} /*end of namespace nde */

#endif

