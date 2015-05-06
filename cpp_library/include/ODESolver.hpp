/**
 * File:   ODESolver.hpp
 * Author: kenobi
 *
 * Created on May 4, 2015, 11:00 PM
 */

#include "Vector.hpp"
#include "Matrix.hpp"

namespace nde {


class ODESolver {
public:
	ODESolver(unsigned int n, double max_error);
	
	Matrix<double> solve(double t0, double tn, double h,
									const Vector<double>& y0) const;

	Vector<double> nextStep(double t, const Vector<double>& y, double h,
								double& err_estimate) const;

	virtual Vector<double> odeSolverDy(double t,
													const Vector<double>& y) const = 0;

private:
	unsigned int _n;
	double _max_error;

	Vector<double> RungeKutta(unsigned int n, unsigned int order, double t,
									const Vector<double>& y, double h) const;

	/* routine with hard-coded Runge Kutta coefficients	*/
	void RungeKuttaCoefficients(unsigned int n, unsigned int order,
		Vector<double>& a, Matrix<double>& b, Vector<double>& c) const;

};


} /*end of namespace nde */

