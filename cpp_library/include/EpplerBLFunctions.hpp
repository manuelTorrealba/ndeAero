/**
 * File:   EpplerBLFunctions.hpp
 * Author: kenobi
 *
 * Created on Jun 16, 2015, 15:00 PM
 */
#ifndef INCLUDE_EPPLERBLFUNCTIONS_HPP
#define INCLUDE_EPPLERBLFUNCTIONS_HPP

namespace nde {

/**
  * Richard Eppler's regression functions from the reference
  * A Computer Program for the Design and Analysis of Low-Speed Airfoils
  * Richard Eppler, Dan M. Somers
  */
class EpplerBLFunctions {
public:
static double calcH12(double h32, bool turbulent);
static double calcCfric(double h32, double Re_d2, bool turbulent);
static double calcCdiss(double h32, double Re_d2, bool turbulent);
static bool turbulentTransition(double h32, double Re_d2, double roughness);
static double calcDeltaCd(double d2, double u, double h12);
static double calcDeltaCl(bool attached, double angle_attack,
								  double length_separation, double slope_te);
};

} /* end of namespace nde */

#endif

