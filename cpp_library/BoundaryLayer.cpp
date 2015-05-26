/**
 * File:   BoundaryLayer.cpp
 * Author: kenobi
 *
 * Created on Apr 25, 2015, 11:00 PM
 */

#include "BoundaryLayer.hpp"
#include <cmath>

namespace nde {

BoundaryLayer::BoundaryLayer(double Re, double U_inf, double l,
				const Interpolator1D& U, const Interpolator1D& dU,
				const Interpolator1D& V0
			) : ODESolver(4, 1.), _Re(Re), _U_inf(U_inf), _l(l),
				_U(U), _dU(dU), _V0(V0) {
	;
}

Vector<double> BoundaryLayer::odeSolverDy(double x,
									const Vector<double>& y) const {

	// outer flow functions
	double u = _U(x);
	double du = _dU(x);
	double v0 = _V0(x);

	Vector<double> dy(2);

	double d2 = y(0);
	double d3 = y(1);
	double H32 = d3 / d2;
	double Re_d2 = _Re * _U(x) * d2 / (_U_inf * _l);

	dy(0) = Cf(H32, Re_d2) + v0 / u - (2.0 + H12(H32, Re_d2)) * du / u * d2;
	dy(1) = Cd(H32, Re_d2) + v0 / u - 3.0 * du / u * d3;

	return dy;

}

double BoundaryLayer::Cf(double H32, double Re_d2) const {

	double eps;

	if (!turbulentTransition(H32, Re_d2)) {
		// laminar boundary layer law

		if (H32 >= 1.57258)
			eps = 1.372391 - 4.226253 * H32 + 2.221687 * H32 * H32;
		else if (H32 >= 1.51509) {
			double h12 = H12(H32, Re_d2);
			eps = 2.512589 - 1.686095 * h12 + 0.391541 * h12 * h12
					  - 0.031729 * h12 * h12 * h12;
		}
		else
			eps = 0.0; // laminar separation!

		return eps / Re_d2;

	} else {
		// turbulent boundary layer laws

		if (H32 >= 1.46) {
			double h12 = H12(H32, Re_d2);
			eps = std::pow((h12 - 1.) * Re_d2, -0.232) * std::exp(-1.26 * h12);	
			return 0.045716 * eps;
		} else
			return 0.0; // turbulent separation!

	}

}

double BoundaryLayer::Cd(double H32, double Re_d2) const {

	if (!turbulentTransition(H32, Re_d2)) {
		// laminar boundary layer law
		if (H32 >= 1.51509) {
			double D = 7.853976 - 10.260551 * H32 + 3.418898 * H32 * H32;
			return 2.0 * D / Re_d2;
		} else
			return 0.0; // laminar separation!
	} else {	// turbulent boundary layer laws
		if (H32 >= 1.46)
			return 0.01 * std::pow((H12(H32, Re_d2) - 1.) * Re_d2 ,-0.1666667);
		else
			return 0.0; // turbulent separation!
	}

}

double BoundaryLayer::H12(double h32, double Re_d2) const {

	if (!turbulentTransition(h32, Re_d2)) {
		// laminar boundary layer law
		if (h32 >= 1.57258)
			return 79.870845 - 89.582142 * h32 + 25.715786 * h32 * h32;
		else if (h32 >= 1.51509)
			return 4.02922 - (583.60182 - 724.55916 * h32
					  + 227.18220 * h32 * h32)
					  * std::sqrt(h32 - 1.51509);
		else
			return 0.0; 	// laminar separation!

	} else { // turbulent boundary layer laws
		if (h32 >= 1.46)
			return (11. * h32 + 15.) / (48. * h32 - 59.);
		else
			return 0.0; // turbulent separation!
	}

}

bool BoundaryLayer::turbulentTransition(double H32, double Re_d2) const {
	return std::log(Re_d2) >= 18.4 * H32 - 21.74 - 0.36 * _roughness;
}

} /* end of namespace nde */

