/**
 * File:   ODESolver.cpp
 * Author: kenobi
 *
 * Created on May 4, 2015, 11:00 PM
 */

#include "ODESolver.hpp"


namespace nde {

ODESolver::ODESolver(unsigned int n, double max_error) :
						_n(n), _max_error(max_error) {
	;
}


Matrix<double> ODESolver::solve(double t0, double tn, unsigned int n_steps,
										const Vector<double>& y0) const {

	double h = (tn-t0)/double(n_steps);
	Matrix<double> y_sol(y0.size(), n_steps + 1);

	// initial condition
	Vector<double> y1 = y0;
	for (unsigned int j = 0; j < y0.size(); ++j)
		y_sol(j,0) = y1(j);

	// evolve in time
	for (unsigned int i = 1; i < n_steps + 1; ++i) {
		double err_estimate;
		double t = t0 + h*(i-1);
		// apply the jump conditions
		y1 = y1 + odeSolverJumpy(t, y1);
		// evolve from t, t+h in time
		Vector<double> y1_next = nextStep(t, y1, h, err_estimate);
		y1 = y1_next;
		for (unsigned int j = 0; j < y0.size(); ++j)
			y_sol(j,i) = y1(j);
	}

	return y_sol;

}


Vector<double> ODESolver::nextStep(double t, const Vector<double>& y,
									  double h, double& err_estimate) const {

	if (_n == 1) { // this is the explicit Euler

		return y + 	odeSolverDy(t, y) * h;

	} else if (_n == 4) { // fourth order integration

		// first, normal evaluation
		Vector<double> y1 = RungeKutta(_n, 4, t, y, h);

		// now, two evaluations with half the step size
		// to get an estimation of error size
		double hh = 0.5 * h;
		Vector<double> yh = RungeKutta(_n, 4, t, y, hh);
		double th = t + hh;
		Vector<double> y1h = RungeKutta(_n, 4, th, yh, hh);

		err_estimate = (y1 - y1h).norm();

		return y1;

	} else if (_n == 6) { // fifth order integration


		Vector<double> y1 = RungeKutta(_n, 4, t, y, h);
		Vector<double> y1h = RungeKutta(_n, 5, t, y, h);

		err_estimate = (y1 - y1h).norm();

		return y1h;

	} else {
		// TODO: throw errortime
	}

}


Vector<double> ODESolver::RungeKutta(unsigned int n, unsigned int order,
					double t, const Vector<double>& y, double h) const {

	Vector<double> a;
	Matrix<double> b;
	Vector<double> c;

	RungeKuttaCoefficients(n, order, a, b, c);

	Matrix<double> k(y.size(), int(n));

	Vector<double> y1 = y; // return vector

	for (unsigned int i = 0; i < n; ++i) {

		Vector<double> ys = y;
		for (unsigned int j = 0; j < i; ++j) {
			for (unsigned int l = 0; l < y.size(); ++l)
				ys(l) = ys(l) + k(l,j) * b(i,j);
		}

		Vector<double>	k_aux = odeSolverDy(t + a(i) * h, ys) * h;

		y1 = y1 + k_aux * c(i);

		for (unsigned int l = 0; l < y.size(); ++l)
			k(l,i) = k_aux(l);

	}

	return y1;

}

void ODESolver::RungeKuttaCoefficients(unsigned int n, unsigned int order,
				Vector<double>& a, Matrix<double>& b, Vector<double>& c) const {

	a.resize(n);
	b.resize(n, n-1);
	b.fill(0.0);
	c.resize(n);

	switch (n) {
	case 4: // typical fourth order RK
	{
		c(0) = 1. / 6.;	a(0) = 0.0;
		c(1) = 1. / 3.;	a(1) = 0.5;
		c(2) = 1. / 3.;	a(2) = 0.5;
		c(3) = 1. / 6.;	a(3) = 1.0;

		b(1,0) = 0.5;
		b(2,1) = 0.5;
		b(3,2) = 1.0;

		break;
	}
	case 6: // six coefficients RK
	{
		a(0) = 0.0;
		a(1) = 0.2;
		a(2) = 0.3;
		a(3) = 0.6;
		a(4) = 1.;
		a(5) = 7. / 8.;

		b(1,0) = 0.2;
		b(2,0) = 3. / 40. ; b(2,1) = 9. / 40.;
		b(3,0) = 0.3; b(3,1) = -0.9; b(3,2) = 6. / 5.;
		b(4,0) = -11. / 54.; b(4,1) = 2.5;
		b(4,2) = -70. / 27. ; b(4,3) = 35. / 27.;
		b(5,0) = 1631. / 55296.; b(5,1) = 175. / 512.; b(5,2) = 575. / 13824.;
		b(5,3) = 44275. / 110592.; b(5,4) = 253. / 4096.;

		if (order == 5) {
			c(0) = 37. / 378.;
			c(1) = 0.;
			c(2) = 250. / 621.;
			c(3) = 125. / 594.;
			c(4) = 0.;
			c(5) = 512. / 1771.;
		} else if (order == 4) {
			c(0) = 2825. / 27648.;
			c(1) = 0.;
			c(2) = 18575. / 48384.;
			c(3) = 13525. / 55296.;
			c(4) = 277. / 14336.;
			c(5) = 0.25;
		} else {
			//TODO: add a runtimeerror.
		}

		break;
	}
	default:
	{
		break;	//TODO: add a runtimeerror.
	}
	}
	
}

} /*end of namespace nde */

