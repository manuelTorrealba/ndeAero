/**
 * File:   AirfoilBoundaryLayer.hpp
 * Author: kenobi
 *
 * Created on Jun 16, 2015, 12:00 PM
 */

#ifndef INCLUDE_AIRFOILBOUNDARYLAYER_HPP
#define INCLUDE_AIRFOILBOUNDARYLAYER_HPP

#include "AerodynamicBody2D.hpp"
#include "BoundaryLayer.hpp"
#include "Vector.hpp"

namespace nde {

/**
  * Application of Boundary Layer Theory to an airfoil with top and bottom
  * surfaces
  */
class AirfoilBoundaryLayer {
public:
	AirfoilBoundaryLayer(double U_inf,
							double chord,
							const AerodynamicBody2D& aerodynamic_body_2D);

	virtual ~AirfoilBoundaryLayer();

	void solveBoundaryLayerFlow();
	double getDeltaCD() const;
	double getDeltaCL() const;

private:
	double _U_inf;
	double _chord;
	double _Re;
	double _angle_attack;
	double _angle_te_top; // top surface trailing edge angle with the x-axis
	double _angle_te_bottom; // bottom surface trailing edge angle with the x-axis
	BoundaryLayer* _boundary_layer_top;
	BoundaryLayer* _boundary_layer_bottom;
	BoundaryLayerSolution* _boundary_layer_sol_top;
	BoundaryLayerSolution* _boundary_layer_sol_bottom;

};

} /* end of namespace nde */

#endif

