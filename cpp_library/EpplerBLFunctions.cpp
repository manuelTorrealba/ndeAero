/**
 * File:   EpplerBLFunctions.cpp
 * Author: kenobi
 *
 * Created on Jun 16, 2015, 15:00 PM
 */
#include "EpplerBLFunctions.hpp"
#include <cmath>
#include <algorithm>

namespace nde {

/******************************************************************************/
/*	Eppler Regression Functions class														*/
/******************************************************************************/
double EpplerBLFunctions::calcH12(double h32, bool turbulent) {

	if (!turbulent) {

		if (h32 >= 1.57258)
			return 79.870845 - 89.58214 * h32 + 25.715784 * h32 * h32;
		else if (h32 >= 1.51509)
			return 4.02922 - (583.60182 - 724.55916 * h32 + 227.18220 * h32 * h32)
								 * std::sqrt(h32 - 1.51509);
		else
			return 0.0; //laminar separation!

	} else

		return (11. * h32 + 15.) / (48. * h32 - 59.);

}


double EpplerBLFunctions::calcCfric(double h32, double Re_d2, bool turbulent) {

	if (!turbulent) {

		double eps;

		if (h32 >= 1.57258)
			eps = 1.372391 - 4.226253 * h32 + 2.221687 * h32 * h32;
		else if (h32 >= 1.51509) {
			double h12 = calcH12(h32, turbulent);
			eps = 2.512589 - 1.686095 * h12 + 0.391541 * h12 * h12
				 - 0.031729 * h12 * h12 * h12;
		} else
			eps = 0.; // laminar separation!

		return eps / Re_d2;

	} else {

		double h12 = calcH12(h32, turbulent);
		return 0.045716 * std::pow((h12 - 1.) * Re_d2 , -0.232)
				* std::exp(-1.260 * h12);

	}
		
		
}

double EpplerBLFunctions::calcCdiss(double h32, double Re_d2, bool turbulent) {

	if (!turbulent) {

		double D = 7.853976 - 10.260551 * h32 + 3.418898 * h32 * h32;
		return 2. * D / Re_d2;

	} else {

		double h12 = calcH12(h32, turbulent);
		return 0.01 * std::pow((h12 - 1.) * Re_d2, -1./6.);

	}

}


bool EpplerBLFunctions::turbulentTransition(double h32, double Re_d2,
														double roughness) {
	return std::log(Re_d2) >= 18.4 * h32 - 21.74 - 0.36 * roughness;
}

double EpplerBLFunctions::calcDeltaCd(double d2, double u, double h12) {
	double hs12 = 0.5 * std::min(h12, 2.5);
	return 2 * d2 * std::pow(u , 2.5 + hs12);
}

double EpplerBLFunctions::calcDeltaCl(bool attached, double angle_attack,
								  double length_separation, double slope_te) {

	if (attached)
		return 0.;
	else
		return - M_PI * length_separation * (slope_te + angle_attack);

}

} /* end of namespace nde */

