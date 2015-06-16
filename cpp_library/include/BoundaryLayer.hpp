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
#include "EpplerBLFunctions.hpp"

namespace nde {

/**
  * Class for "easy" handling Boundary Layer Result Vectors
  */
class BoundaryLayerSolution {
public:
	BoundaryLayerSolution(const Matrix<double>& boundary_layer_ode_solution);
	BoundaryLayerSolution(const BoundaryLayerSolution& boundary_layer_solution);

	virtual ~BoundaryLayerSolution();

	bool getIsAttached() const;
	double getSeparationX() const;
	double calcD1(double x) const;
	double calcD2(double x) const;
	double calcD3(double x) const;
	double calcH12(double x) const;
	double calcH32(double x) const;

private:
	bool _is_attached;
	double _separation_x;
	Interpolator1D* _d1;
	Interpolator1D* _d2;
	Interpolator1D* _d3;
	Interpolator1D* _h12;
	Interpolator1D* _h32;
};

/**
  * Boundary Layer equations from reference
  * A Computer Program for the Design and Analysis of Low-Speed Airfoils
  * Richard Eppler, Dan M. Somers
  */
class BoundaryLayer : public ODESolver {
public:
	BoundaryLayer(double Re, const Interpolator1D& U_ext);

	// main solution function
	BoundaryLayerSolution solve() const;

	// inherited from ODESolver
	virtual Vector<double> odeSolverDy(double x,
												const Vector<double>& y) const;
	virtual Vector<double> odeSolverJumpy(double x,
													const Vector<double>& y) const;

	// get methods
	const Interpolator1D& getUExt() const;

private:
	double _Re;
	double _roughness; // TODO: initialize roughness!
	Interpolator1D _U_ext;
	Vector<double> setBlausiusInitialCondition(double x_0) const;

};

} /* end of namespace nde */

#endif

