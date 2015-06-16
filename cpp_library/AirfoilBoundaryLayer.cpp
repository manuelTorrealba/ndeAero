/**
 * File:   AirfoilBoundaryLayer.cpp
 * Author: kenobi
 *
 * Created on Jun 16, 2015, 12:00 PM
 */

#include "AirfoilBoundaryLayer.hpp"
#include "Constants.hpp"

namespace nde {

/******************************************************************************/
/*	Class AirfoilBoundaryLayer																	*/
/******************************************************************************/

AirfoilBoundaryLayer::AirfoilBoundaryLayer(double U_inf, double chord,
							const AerodynamicBody2D& aerodynamic_body_2D)
						 : _U_inf(U_inf), _chord(chord) {

	// Reynolds number
	_Re = _U_inf * _chord / AIR_KINEMATIC_VISCOSITY;

	// this is the arc-length vector measured from the
	// botton-surface trailing edge and the corresponding air-flow speed.
	Vector<double> d_x = aerodynamic_body_2D.getControlPointsDistances();
	Vector<double> v_x = aerodynamic_body_2D.getVelocities();


	// angle of attack	
	_angle_attack = aerodynamic_body_2D.getAngleAttack();

	// airfoil geometric angles at trailing edge
	// TODO: this has to be done properly with atan2
	_angle_te_top = -std::acos(std::abs(
								aerodynamic_body_2D.getTangents()(d_x.size() - 1)(0)));
	_angle_te_bottom = std::acos(std::abs(
								aerodynamic_body_2D.getTangents()(0)(0)));

	// find stagnation point
	// the velocity function will show a minimum around the stagnation point
	size_t index_smallest_v = 0;
	for (size_t i = 1; i < v_x.size()-1; ++i) {
		if (v_x(i-1) >= v_x(i) && v_x(i+1) >= v_x(i)) {
			index_smallest_v = i;
			break;
		}
	}

	double d_x_smallest_v;
	if (index_smallest_v > 0 && index_smallest_v < v_x.size() - 1) {
		if (v_x(index_smallest_v - 1) < v_x(index_smallest_v + 1)) {
			d_x_smallest_v = 0.5 * (d_x(index_smallest_v - 1)
										 + d_x(index_smallest_v));
			index_smallest_v = index_smallest_v - 1;
		} else
			d_x_smallest_v = 0.5 * (d_x(index_smallest_v + 1)
										 + d_x(index_smallest_v));
	} else
		d_x_smallest_v = d_x(index_smallest_v);

	// arrange velocity vectors for top and bottom surfaces
	Vector<double> x_bottom(index_smallest_v + 2);
	Vector<double> v_bottom(index_smallest_v + 2);
	x_bottom(0) = 0.; // stagnation point
	v_bottom(0) = 0.;
	for (size_t i = 0; i <= index_smallest_v; ++i) {
		x_bottom(i + 1) = d_x_smallest_v - d_x(index_smallest_v - i);
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

	// build external flow curves for top and bottom surfaces
	Interpolator1D U_ext_bottom(x_bottom, v_bottom, SPLINE_MONOTONE);
	Interpolator1D U_ext_top(x_top, v_top, SPLINE_MONOTONE);

	// initialize top and bottom boundary layer objects
	_boundary_layer_bottom = new BoundaryLayer(_Re, U_ext_bottom);
	_boundary_layer_top = new BoundaryLayer(_Re, U_ext_top);

/*
	std::cout << "bottom" << std::endl;
	for (size_t i = 0; i < x_bottom.size(); ++i) {
		std::cout << x_bottom(i) << "," << v_bottom(i) << std::endl;
	}

	std::cout << "top" << std::endl;
	for (size_t i = 0; i < x_top.size(); ++i) {
		std::cout << x_top(i) << "," << v_top(i) << std::endl;
	}
*/

}

//
// destructor
//
AirfoilBoundaryLayer::~AirfoilBoundaryLayer() {

	delete _boundary_layer_bottom;
	delete _boundary_layer_top;
	if (_boundary_layer_sol_top)
		delete _boundary_layer_sol_top;
	if (_boundary_layer_sol_bottom)
		delete _boundary_layer_sol_bottom;

}

void AirfoilBoundaryLayer::solveBoundaryLayerFlow() {

	_boundary_layer_sol_top = new BoundaryLayerSolution
												(_boundary_layer_top->solve());

	_boundary_layer_sol_bottom = new BoundaryLayerSolution
												(_boundary_layer_bottom->solve());

}

double AirfoilBoundaryLayer::getDeltaCD() const {

	// top surface
	double x_top;
	if (_boundary_layer_sol_top->getIsAttached())
		x_top = (_boundary_layer_top->getUExt()).getInterpolationRange()(1);
	else
		x_top = _boundary_layer_sol_top->getSeparationX();

	double u_top = (_boundary_layer_top->getUExt())(x_top);
	double d2_top = _boundary_layer_sol_top->calcD2(x_top);
	double h12_top = _boundary_layer_sol_top->calcH12(x_top);
	double delta_cd_top = EpplerBLFunctions::calcDeltaCd(d2_top, u_top, h12_top);

	// bottom surface
	double x_bottom;
	if (_boundary_layer_sol_bottom->getIsAttached())
		x_bottom = (_boundary_layer_bottom->getUExt()).getInterpolationRange()(1);
	else
		x_bottom = _boundary_layer_sol_bottom->getSeparationX();

	double u_bottom = (_boundary_layer_bottom->getUExt())(x_bottom);
	double d2_bottom = _boundary_layer_sol_bottom->calcD2(x_bottom);
	double h12_bottom = _boundary_layer_sol_bottom->calcH12(x_bottom);
	double delta_cd_bottom = EpplerBLFunctions::calcDeltaCd(d2_bottom, u_bottom,
															h12_bottom);
	return delta_cd_top + delta_cd_bottom;

}


double AirfoilBoundaryLayer::getDeltaCL() const {

	// top surface
	double delta_cl_top;

	if (_boundary_layer_sol_top->getIsAttached())

		delta_cl_top = 0.0;

	else {

		double length_separation_top =
			(_boundary_layer_top->getUExt()).getInterpolationRange()(1)
		 - _boundary_layer_sol_top->getSeparationX();

		delta_cl_top = EpplerBLFunctions::calcDeltaCl(false, _angle_attack,
								  length_separation_top, _angle_te_top);

	}

	// bottom surface
	double delta_cl_bottom;

	if (_boundary_layer_sol_bottom->getIsAttached())

		delta_cl_bottom = 0.0;

	else {

		double length_separation_bottom =
			(_boundary_layer_bottom->getUExt()).getInterpolationRange()(1)
		 - _boundary_layer_sol_bottom->getSeparationX();

		delta_cl_bottom = EpplerBLFunctions::calcDeltaCl(false, _angle_attack,
								  length_separation_bottom, _angle_te_bottom);

	}

	return std::min(delta_cl_top, 0.0) + std::max(delta_cl_bottom, 0.0);

}

} /* end of namespace nde */

