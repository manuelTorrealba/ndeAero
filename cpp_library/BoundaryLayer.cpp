/**
 * File:   BoundaryLayer.cpp
 * Author: kenobi
 *
 * Created on Apr 25, 2015, 11:00 PM
 */

#include "BoundaryLayer.hpp"
#include "Constants.hpp"
#include <cmath>

namespace nde {

BoundaryLayer::BoundaryLayer(double Re, const Interpolator1D& U)
									: ODESolver(1, 1.), _Re(Re), _U(U) {
	_roughness = 0.;
}


Vector<double> BoundaryLayer::odeSolverJumpy(double x,
														const Vector<double>& y) const {

	// outer flow functions
	double u = _U(x);
	double d2 = y(0);
	double d3 = y(1);
	double turbulent_flag = y(2);
	double separation_flag = y(3);

	// check for turbulent regime
	double Re_d2 = _Re * u * d2;
	double h32 = d3 / d2;
	bool is_turbulent = (turbulent_flag < -0.5 || turbulentTransition(h32, Re_d2));

	// check for the flow still attached
	bool is_attached = !(separation_flag < -0.5);

	// if it is still attached check for conditions
	// the flow being attached or not at this stage
	if (is_attached) {

		// Continue with turbulent equations.
		// In practise what probably will happen is that
		// flow will be reattached at a later section.
		if (!is_turbulent && h32 < 1.51509) is_turbulent = true; 

		//separation surely happens.
		if (is_turbulent && h32 < 1.46) is_attached = false;

	}

	// apply a negative jump for the turbulent and attached variables.
	Vector<double> Jump_y(4, 0.0);
	if (is_turbulent) Jump_y(2) = -1.;
	if (!is_attached) Jump_y(3) = -1.;

	return Jump_y;

}

Vector<double> BoundaryLayer::odeSolverDy(double x,
									const Vector<double>& y) const {

	bool is_turbulent = (y(2) < -0.5);
	bool is_attached = !(y(3) < -0.5);

	Vector<double> dy(4, 0.0);
	if (is_attached) {

		// outer flow functions
		double u = _U(x);
		double du = _U(1, x);

		// boundary layer current status
		double d2 = y(0);
		double d3 = y(1);
		double Re_d2 = _Re * u * d2;
		double h32 = d3 / d2;

		double h12 = calcH12(h32, is_turbulent);
		double cf = calcCf(h32, Re_d2, is_turbulent);
		double cd = calcCd(h32, Re_d2, is_turbulent);

		dy(0) = cf - (2.0 + h12) * du / u * d2;
		dy(1) = cd - 3.0 * du / u * d3;

	}

	return dy;

}


double BoundaryLayer::calcH12(double h32, bool turbulent) const {

	if (!turbulent) {

		if (h32 >= 1.57258)
			return 79.870845 - 89.58214 * h32 + 25.715784 * h32 * h32;
		else if (h32 >= 1.51509)
			return 4.02922 - (583.60182 - 724.55916 * h32 + 227.18220 * h32 * h32)
								 * std::sqrt(h32 - 1.51509);
		else
			return 0.0; //laminar separation!

	} else {

		return (11. * h32 + 15.) / (48. * h32 - 59.);

	}

}


double BoundaryLayer::calcCf(double h32, double Re_d2, bool turbulent) const {

	if (!turbulent) {

		double eps;

		if (h32 >= 1.57258)
			eps = 1.372391 - 4.226253 * h32 + 2.221687 * h32 * h32;
		else if (h32 >= 1.51509) {
			double h12 = calcH12(h32, turbulent);
			eps = 2.512589 - 1.686095 * h12 + 0.391541 * h12 * h12
				 - 0.031729 * h12 * h12 * h12;
		} else {
			eps = 0.; // laminar separation!
		}

		return eps / Re_d2;

	} else {

		double h12 = calcH12(h32, turbulent);
		return 0.045716 * std::pow((h12 - 1.) * Re_d2 , -0.232)
				* std::exp(-1.260 * h12);

	}
		
		
}


double BoundaryLayer::calcCd(double h32, double Re_d2, bool turbulent) const {

	if (!turbulent) {

		double D = 7.853976 - 10.260551 * h32 + 3.418898 * h32 * h32;
		return 2. * D / Re_d2;

	} else {

		double h12 = calcH12(h32, turbulent);
		return 0.01 * std::pow((h12 - 1.) * Re_d2, -1./6.);

	}

}


bool BoundaryLayer::turbulentTransition(double h32, double Re_d2) const {
	return std::log(Re_d2) >= 18.4 * h32 - 21.74 - 0.36 * _roughness;
}


/******************************************************************************/
/*	Class AirfoilBoundaryLayer																	*/
/******************************************************************************/

AirfoilBoundaryLayer::AirfoilBoundaryLayer(double U_inf, double chord,
							const Vector<double>& d_x, // this is the arc-length vector
																// measured from the
																// botton-surface trailing edge
							const Vector<double>& v_x) : _U_inf(U_inf), _chord(chord) {

	// Reynolds number
	_Re = _U_inf * _chord / AIR_KINEMATIC_VISCOSITY;

	std::cout << "Reynolds number = " << _Re << std::endl;

	// find stagnation point
	size_t index_smallest_v = v_x.smallestAbsIndex();
	double d_x_smallest_v;
	if (index_smallest_v > 0 && index_smallest_v < v_x.size() - 1) {
		if (v_x(index_smallest_v - 1) > v_x(index_smallest_v))
			d_x_smallest_v = 0.5 * (d_x(index_smallest_v)
										 + d_x(index_smallest_v + 1));
		else
			d_x_smallest_v = 0.5 * (d_x(index_smallest_v - 1)
										 + d_x(index_smallest_v));
	} else
		d_x_smallest_v = d_x(index_smallest_v);

	// arrange velocity vectors for top and bottom surfaces
	Vector<double> x_bottom(index_smallest_v + 2);
	Vector<double> v_bottom(index_smallest_v + 2);
	x_bottom(0) = 0.; // stagnation point
	v_bottom(0) = 0.;
	for (size_t i = 0; i <= index_smallest_v; ++i) {
		x_bottom(i + 1) = d_x(index_smallest_v - i);
		v_bottom(i + 1) = v_x(index_smallest_v - i);
	}

	Vector<double> x_top(v_x.size() - index_smallest_v);
	Vector<double> v_top(v_x.size() - index_smallest_v);
	x_top(0) = 0.; // stagnation point
	v_top(0) = 0.;
	for (size_t i = index_smallest_v + 1; i < v_x.size(); ++i) {
		x_top(i - index_smallest_v) = d_x(i) - d_x_smallest_v;
		v_top(i - index_smallest_v) = v_x(i);
	}


	Interpolator1D U_bottom(x_bottom, v_bottom, SPLINE_MONOTONE);
	Interpolator1D U_top(x_top, v_top, SPLINE_MONOTONE);
//	Interpolator1D U_top(Vector<double>(1,1.0), Vector<double>(1,1.0),
//								SPLINE_MONOTONE);

	// initialize top and bottom boundary layer objects
	_boundary_layer_bottom = new BoundaryLayer(_Re, U_bottom);
	_boundary_layer_top = new BoundaryLayer(_Re, U_top);

	// solve PDEs numerically
	/*
	Vector<double> y0_bottom(2);
	y0_bottom(0) = 0.29004 / std::sqrt(_Re * U_bottom(1,0.0)); // blausius solution
	y0_bottom(1) = 1.61998 * y0_bottom(0);

	Matrix<double> y_bottom = _boundary_layer_bottom->solve
										(0., x_bottom.last_v(), x_bottom.last_v() / 10.,
										y0_bottom);
	*/

	// start with blausius solution
	double x_0_top = 0.01;
	double Re_0_top = _Re * U_top(x_0_top) * x_0_top;
	Vector<double> y0_top(4);
	y0_top(0) = 0.29004 * x_0_top / std::sqrt(Re_0_top); // blausius solution
	y0_top(1) = 1.61998 * y0_top(0); // blausius solution
	y0_top(2) = 0.; // assume flow initially laminar
	y0_top(3) = 0.; // assume flow initially attached

	Matrix<double> y_top = _boundary_layer_top->solve
									(x_0_top, x_top.last_v(), std::pow(2, 8), y0_top);

	std::cout<< "boundary layer top" << std::endl;
	for (size_t i = 0; i < 21; ++i) {
		double s = i * x_top.last_v() / 20.;
		std::cout << s << ", " << y_top(0, i) << ", " << y_top(1, i) <<
								", " << y_top(2, i) << ", " << y_top(3, i) << std::endl;
	}


}

AirfoilBoundaryLayer::~AirfoilBoundaryLayer() {

	delete _boundary_layer_bottom;
	delete _boundary_layer_top;

}


} /* end of namespace nde */

