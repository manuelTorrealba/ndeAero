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

class BoundaryLayer : public ODESolver {
public:
BoundaryLayer(double Re, double U_inf, double l,
				const Interpolator1D& U, const Interpolator1D& dU,
				const Interpolator1D& V0);

virtual Vector<double> odeSolverDy(double x, const Vector<double>& y) const;

private:
double _Re;
double _U_inf;
double _l;
double _roughness;
Interpolator1D _U;
Interpolator1D _dU;
Interpolator1D _V0;

double Cf(double H32, double Re_d2) const;
double Cd(double H32, double Re_d2) const;
double H12(double H32, double Re_d2) const;
bool turbulentTransition(double H32, double Re_d2) const;
};

} /* end of namespace nde */

#endif

