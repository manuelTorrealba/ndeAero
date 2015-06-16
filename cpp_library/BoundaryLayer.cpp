/**
 * File:   BoundaryLayer.cpp
 * Author: kenobi
 *
 * Created on Apr 25, 2015, 11:00 PM
 */

#include "BoundaryLayer.hpp"
#include <cmath>

namespace nde {

BoundaryLayer::BoundaryLayer(double Re, const Interpolator1D& U_ext)
									: ODESolver(1, 1.), _Re(Re), _U_ext(U_ext) {
	_roughness = 0.;
}

BoundaryLayerSolution BoundaryLayer::solve() const {

	//
	// ODE solution for the boundary layer equations
	//

	// high number of steps -> solution needs small increments for ODE solver
	//									stability
	double n_steps = std::pow(2, 14);

	// Initial solution for ODE solver:
	// start with blausius solution for a small x_0 = 0.01
	double x_0 = 0.01;

	return BoundaryLayerSolution(this->ODESolver::solve(
														x_0,
														_U_ext.getInterpolationRange()(1),
														n_steps,
														setBlausiusInitialCondition(x_0)));

}

Vector<double> BoundaryLayer::odeSolverJumpy(double x,
														const Vector<double>& y) const {

	// this routine applies a negative jump for
	// the turbulent (component 2 of the result vector)
	// and attached (component 3 of the return vector) signal variables.
	Vector<double> Jump_y(4, 0.0);

	// check for the flow still attached
	double separation_flag = y(3);
	bool is_attached = !(separation_flag < -0.5);

	// if it is still attached check for conditions
	// the flow being attached or not at this stage
	if (is_attached) {

		// outer flow functions
		double u = _U_ext(x);
		double d2 = y(0);
		double d3 = y(1);
		double turbulent_flag = y(2);

		// check for turbulent regime

		bool is_turbulent = (turbulent_flag < -0.5);

		// check transition from laminar to turbulent, as a function
		// of Re_d2 and h32.
		double Re_d2 = _Re * u * d2;
		double h32 = d3 / d2;

		// In case of laminar flow separation (h32 < 1.51509),
		// continue with turbulent boundary layer equations.
		// In practise what probably will happen is that
		// flow will be reattached at a later section.
		if (!is_turbulent &&
			 (EpplerBLFunctions::turbulentTransition(h32, Re_d2, _roughness)
			  || h32 < 1.51509)) {

			is_turbulent = true;
			// jump by -1 to indicate turbulent transition
			// at this point
			Jump_y(2) = -1.;

		}

		//separation surely happens.
		if (is_turbulent && h32 < 1.46) {
			is_attached = false;
			// jump by -1 to indicate flow separation
			// at this point
			Jump_y(3) = -1.;
		}

	}

	return Jump_y;

}

Vector<double> BoundaryLayer::odeSolverDy(double x,
													const Vector<double>& y) const {

	bool is_turbulent = (y(2) < -0.5);
	bool is_attached = !(y(3) < -0.5);

	Vector<double> dy(4, 0.0);
	if (is_attached) {

		// outer flow functions
		double u = _U_ext(x);
		double du = _U_ext(1, x);

		// boundary layer current status
		double d2 = y(0);
		double d3 = y(1);
		double Re_d2 = _Re * u * d2;
		double h32 = d3 / d2;

		double h12 = EpplerBLFunctions::calcH12(h32, is_turbulent);
		double cf = EpplerBLFunctions::calcCfric(h32, Re_d2, is_turbulent);
		double cd = EpplerBLFunctions::calcCdiss(h32, Re_d2, is_turbulent);

		dy(0) = cf - (2.0 + h12) * du / u * d2;
		dy(1) = cd - 3.0 * du / u * d3;

	}

	return dy;

}


Vector<double> BoundaryLayer::setBlausiusInitialCondition(double x_0) const {

	double u_0 = _U_ext(x_0);
	double Re_0 = _Re * u_0 * x_0;
	Vector<double> y_0(4);
	// blausius solution for d2
	y_0(0) = 0.29004 * x_0 / std::sqrt(Re_0);
	// blausius solution for d3
	y_0(1) = 1.61998 * y_0(0);
	y_0(2) = 0.; // assume flow initially laminar
	y_0(3) = 0.; // assume flow initially attached

	return y_0;

}

const Interpolator1D& BoundaryLayer::getUExt() const {
	return _U_ext;
}

/******************************************************************************/
/* class BoundaryLayerSolution																*/
/******************************************************************************/
BoundaryLayerSolution::BoundaryLayerSolution
						(const Matrix<double>& boundary_layer_ode_solution) {

	size_t num_sol_points = boundary_layer_ode_solution.numcols();

	Vector<double> x(num_sol_points);
	Vector<double> d1(num_sol_points);
	Vector<double> d2(num_sol_points);
	Vector<double> d3(num_sol_points);
	Vector<double> h12(num_sol_points);
	Vector<double> h32(num_sol_points);

	for (unsigned int i = 0; i < num_sol_points; ++i) {
		x(i) = boundary_layer_ode_solution(0, i);
		d2(i) = boundary_layer_ode_solution(1, i);
		d3(i) = boundary_layer_ode_solution(2, i);
		h32(i) = d3(i) / d2(i);
		h12(i) = EpplerBLFunctions::calcH12(h32(i),
												boundary_layer_ode_solution(3, i) < -0.5);
		d1(i) = h12(i) * d2(i);
	}

	_d1 = new Interpolator1D(x, d1, SPLINE_MONOTONE);
	_d2 = new Interpolator1D(x, d2, SPLINE_MONOTONE);
	_d3 = new Interpolator1D(x, d3, SPLINE_MONOTONE);
	_h12 = new Interpolator1D(x, h12, SPLINE_MONOTONE);
	_h32 = new Interpolator1D(x, h32, SPLINE_MONOTONE);

	_is_attached = !(boundary_layer_ode_solution(4, num_sol_points - 1) < -0.5);

	if (!_is_attached) {

		// look for the separation point
		for (size_t i = 1; i < num_sol_points; ++i) {

			// when separation is found, determine separation
			// coordinate by imposing h32 = 1.46 and linear interpolation
			if (boundary_layer_ode_solution(4, i) < -0.5) {
				double x0 = x(i - 2);
				double x1 = x(i - 1);
				double y0 = h32(i - 2);
				double y1 = h32(i - 1);
//				std::cout << i << "," << x0 << "," << x1 << "," << y0 << "," << y1 << std::endl;
				_separation_x = x0 + (x1 - x0) / (y1 - y0) * (1.46 - y0);
				break;
			}

		}

	}

}

BoundaryLayerSolution::BoundaryLayerSolution
						(const BoundaryLayerSolution& boundary_layer_solution) {

	_is_attached = boundary_layer_solution._is_attached;
	_separation_x = boundary_layer_solution._separation_x;
	_d1 = new Interpolator1D(*boundary_layer_solution._d1);
	_d2 = new Interpolator1D(*boundary_layer_solution._d2);
	_d3 = new Interpolator1D(*boundary_layer_solution._d3);
	_h12 = new Interpolator1D(*boundary_layer_solution._h12);
	_h32 = new Interpolator1D(*boundary_layer_solution._h32);

}

BoundaryLayerSolution::~BoundaryLayerSolution() {

	delete _d1;
	delete _d2;
	delete _d3;
	delete _h12;
	delete _h32;

}

bool BoundaryLayerSolution::getIsAttached() const {
	return _is_attached;
}

double BoundaryLayerSolution::getSeparationX() const {
	return _separation_x;
}

double BoundaryLayerSolution::calcD1(double x) const {
	return (*_d1)(x);
}

double BoundaryLayerSolution::calcD2(double x) const {
	return (*_d2)(x);
}

double BoundaryLayerSolution::calcD3(double x) const {
	return (*_d3)(x);
}

double BoundaryLayerSolution::calcH12(double x) const {
	return (*_h12)(x);
}

double BoundaryLayerSolution::calcH32(double x) const {
	return (*_h32)(x);
}

} /* end of namespace nde */

