/*
 * File:   Wing.hpp
 * Author: kenobi
 *
 * Created on August 19, 2014, 14:37 PM
 */

#ifndef WING_HPP
#define WING_HPP

#include"Airfoil.hpp"
#include"Panel.hpp"

namespace nde {

class Wing {
public:
	//constructor for a trapezoidal wing
	Wing(const Airfoil& airfoil_root_in, const Airfoil& airfoil_tip_in,
			double wing_span_in, double dihedral_angle_in, // in radians
			double twist_angle_in, //in radians
			double swept_angle_in //in radians
			);

	Vector<Panel3D> getPanelsSimpleSurface(int leftright, int updown, 
					double density_x, double density_y) const;

	Vector<Panel3D> getPanels(double density_x, double density_y) const;

	double getSurface() const;
	double getMeanChord() const;
	Vector<double> getPoint(double x_unit, double y_unit, int top_bottom) const;

private:
	Airfoil airfoil_root;
	Airfoil airfoil_tip;
	double wing_span;
	double dihedral_angle;
	double twist_angle;
	double swept_angle;

	Vector<double> getPointOnWingPlane(double x_unit, double y_unit) const;
	double getAirfoilTop(double x_unit, double y_unit) const;
	double getAirfoilBottom(double x_unit, double y_unit) const;

};

} //end of namespace nde

#endif/*WING_HPP*/
