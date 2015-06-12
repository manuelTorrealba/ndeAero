/**
 * File:   BoundaryLayer.hpp
 * Author: kenobi
 *
 * Created on May 7, 2015, 15:00 PM
 */
#ifndef INCLUDE_BOUNDARYLAYER_HPP
#define INCLUDE_BOUNDARYLAYER_HPP

#include "ODESolver.hpp"
#include "Interpolator.hpp"

namespace nde {

/**
 *
 * Boundary Layer equations from reference
 * A Computer Program for the Design and Analysis of Low-Speed Airfoils
 * Richard Eppler, Dan M. Somers
 */

class BoundaryLayer : public ODESolver {
public:
	BoundaryLayer(double Re, const Interpolator1D& U);

	virtual Vector<double> odeSolverDy(double x,
												const Vector<double>& y) const;
	virtual Vector<double> odeSolverJumpy(double x,
													const Vector<double>& y) const;

private:
	double _Re;
	double _roughness; // TODO: initialize roughness!
	Interpolator1D _U;

	double calcH12(double h32, bool turbulent) const;
	double calcCf(double h32, double Re_d2, bool turbulent) const;
	double calcCd(double h32, double Re_d2, bool turbulent) const;

	bool turbulentTransition(double h, double Re_theta) const;

};



class AirfoilBoundaryLayer {
public:
	AirfoilBoundaryLayer(double U_inf,
							double chord,
							const Vector<double>& d_x, // this is the arc-length vector
																// measured from the
																// botton-surface trailing edge
							const Vector<double>& v_x);

	virtual ~AirfoilBoundaryLayer();

	void solveBoundaryLayerFlow() const;
	double getDeltaCD() const;
	double getDeltaCL() const;

private:
	double _U_inf;
	double _chord;
	double _Re;
	BoundaryLayer* _boundary_layer_top;
	BoundaryLayer* _boundary_layer_bottom;
};

} /* end of namespace nde */

#endif

