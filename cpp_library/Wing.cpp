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
	Vector<Panel3D> z(num_panels_chord * num_panels_span * 4);

	double dx_root = airfoil_root.getChord() / double(num_panels_chord);
	double dx_tip = airfoil_tip.getChord() / double(num_panels_chord);
	double dy_span = 0.5 * wing_span / double(num_panels_span);

	for (int i = 0; i < num_panels_chord; ++i) {

		double x_root_i = double(i) * dx_root;
		double x_root_i1 = double(i + 1) * dx_root;

		double x_tip_i = double(i) * dx_tip;
		double x_tip_i1 = double(i + 1) * dx_tip;

		double z_root_top_i = airfoil_root.getPointTop(x_root_i);
		double z_root_top_i1 = airfoil_root.getPointTop(x_root_i1);

		double z_tip_top_i = airfoil_tip.getPointTop(x_root_i);
		double z_tip_top_i1 = airfoil_tip.getPointTop(x_root_i1);

		double z_root_bottom_i = airfoil_root.getPointBottom(x_root_i);
		double z_root_bottom_i1 = airfoil_root.getPointBottom(x_root_i1);

		double z_tip_bottom_i = airfoil_tip.getPointBottom(x_root_i);
		double z_tip_bottom_i1 = airfoil_tip.getPointBottom(x_root_i1);

		for (int j = 0; j < num_panels_span; ++j) {

			double alpha = double(j) / double(num_panels_span);
			double alpha1 = double(j + 1) / double(num_panels_span);

			Vector<double> x1(3);
			x1(0) = x_root_i * alpha + x_tip_i * (1.0 - alpha);
			x1(1) = 0.5 * wing_span * alpha;
			x1(2) = z_root_top_i * alpha + z_tip_top_i * (1.0 - alpha);

			Vector<double> x2(3);
			x2(0) = x_root_i1 * alpha + x_tip_i1 * (1.0 - alpha);
			x2(1) = 0.5 * wing_span * alpha;
			x2(2) = z_root_top_i1 * alpha + z_tip_top_i1 * (1.0 - alpha);

			Vector<double> x3(3);
			x3(0) = x_root_i1 * alpha1 + x_tip_i1 * (1.0 - alpha1);
			x3(1) = 0.5 * wing_span * alpha1;
			x3(2) = z_root_top_i1 * alpha1 + z_tip_top_i1 * (1.0 - alpha1);

			Vector<double> x4(3);
			x4(0) = x_root_i * alpha1 + x_tip_i * (1.0 - alpha1);
			x4(1) = 0.5 * wing_span * alpha1;
			x4(2) = z_root_top_i * alpha1 + z_tip_top_i * (1.0 - alpha1);

			z(i + j * num_panels_chord).setPoints(x1, x2, x3, x4);

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

}/*end namespace nde */
