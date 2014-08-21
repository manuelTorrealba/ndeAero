/*
 * File:   Wing.cpp
 * Author: kenobi
 *
 * Created on August 19, 2014, 14:37 PM
 */

#include"Wing.hpp"

#include<cmath>

namespace nde {

Wing::Wing(const Airfoil& airfoil_root_in, const Airfoil& airfoil_tip_in,
		double wing_span_in, double dihedral_angle_in, double twist_angle_in,
		double swept_angle_in) :
		airfoil_root(airfoil_root_in), airfoil_tip(airfoil_tip_in), wing_span(
				wing_span_in), dihedral_angle(dihedral_angle_in), twist_angle(
				twist_angle_in), swept_angle(swept_angle_in) {
	;
}

Vector<Panel3D> Wing::getPanels(double density_x, double density_y) const {

	int num_panels_chord = std::floor(1.0 / density_x);
	int num_panels_span = std::floor(1.0 / density_y);
	double dx = 1.0 / double(num_panels_chord);
	double dy = 1.0 / double(num_panels_span);

	int N = num_panels_chord * num_panels_span;
	Vector<Panel3D> z(4 * N);

	for (int i = 0; i < num_panels_chord; ++i) {

		double xi = dx * double(i) / double(num_panels_chord);
		double xi1 = dx * double(i + 1) / double(num_panels_chord);

		for (int j = 0; j < num_panels_span; ++j) {

			double yi = dy * double(j) / double(num_panels_span);
			double yi1 = dy * double(j + 1) / double(num_panels_span);

			// panel in the left wing, top surface
			z(i + j * num_panels_chord).setPoints(getPoint(xi, yi, 1),
					getPoint(xi1, yi, 1), getPoint(xi1, yi1, 1),
					getPoint(xi, yi1, 1));

			// panel in the left wing bottom surface
			z(N + i + j * num_panels_chord).setPoints(getPoint(xi, yi, -1),
					getPoint(xi, yi1, -1), getPoint(xi1, yi1, -1),
					getPoint(xi1, yi, -1));

			// panel in the right wing, top surface
			z(2 * N + i + j * num_panels_chord).setPoints(getPoint(xi, -yi, 1),
					getPoint(xi, -yi1, 1), getPoint(xi1, -yi1, 1),
					getPoint(xi1, -yi, 1));

			// panel in the right wing bottom surface
			z(3 * N + i + j * num_panels_chord).setPoints(getPoint(xi, -yi, -1),
					getPoint(xi1, -yi, -1), getPoint(xi1, -yi1, -1),
					getPoint(xi, -yi1, -1));

		}

	}

	return z;

}

double Wing::getSurface() const {
	return getMeanChord() * wing_span;
}

double Wing::getMeanChord() const {
	return 0.5 * (airfoil_tip.getChord() + airfoil_root.getChord());
}

Vector<double> Wing::getPoint(double x_unit, double y_unit,
		int top_bottom) const {

	Vector<double> x = getPointOnWingPlane(x_unit, y_unit);

	double z = 0.0;

	switch (top_bottom) {
	case 1: {
		z = getAirfoilTop(x_unit, y_unit);
		break;
	}
	case -1: {
		z = getAirfoilBottom(x_unit, y_unit);
		break;
	}
	default:
		break;
	}

	x(2) = x(2) + z * top_bottom;

	return x;

}

// this routine maps [0,1]x[0,1] rectangle on the Wing's surface place.
Vector<double> Wing::getPointOnWingPlane(double x_unit, double y_unit) const {

	double y_sign = 1.0;
	if (y_unit < 0.0)
		y_sign = -1.0;

	y_unit = std::abs(y_unit);

	// chord for this section
	double chord = airfoil_root.getChord() * (1.0 - y_unit)
			- airfoil_tip.getChord() * y_unit;

	Vector<double> x(3);
	// x-coordinate
	x(0) = y_unit * 0.5 * wing_span * std::sin(swept_angle)
			+ (x_unit - 0.25) * chord;
	// y-coordinate
	x(1) = y_sign * (y_unit * 0.5 * wing_span * std::cos(swept_angle));
	// z-coordinate
	x(2) = y_unit * std::sin(dihedral_angle)
			+ (x_unit - 0.25) * y_unit * twist_angle;

	return x;

}

double Wing::getAirfoilTop(double x_unit, double y_unit) const {

	return airfoil_root.getPointTop(x_unit) * (1.0 - y_unit)
			+ airfoil_tip.getPointTop(x_unit) * y_unit;

}

double Wing::getAirfoilBottom(double x_unit, double y_unit) const {

	return airfoil_root.getPointBottom(x_unit) * (1.0 - y_unit)
			+ airfoil_tip.getPointBottom(x_unit) * y_unit;

}

}/*end namespace nde */
