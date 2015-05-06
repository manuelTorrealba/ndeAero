/**
 * File:   BoundaryLayer.cpp
 * Author: kenobi
 *
 * Created on Apr 25, 2015, 11:00 PM
 */

#include "BoundaryLayer.hpp"
#include <cmath>

namespace nde {


	Vector<double> BoundaryLayer::odeSolverDy(double x) const {

		double d2 = D2(x);
		double d3 = D3(x);
		double h32 = d3 / d2;
		double h12 = H12(h32);
		double u = U(x);
		double du = dU(x);
		double v0 = V0(x);

		double dd2 = Cf(x, Eps(h32), d2) + v0 / u
						 - (2.0 + h12) * du / u * d2;
		double dd3 = Cd(x, D(h32), d2) + v0 / u
						 - 3.0 * du / u * d3;

	}

	double BoundaryLayer::Cf(double x, double eps, double d2) const {
		return eps / (_Re * U(x) / _U_inf * d2);
	}

	double BoundaryLayer::Cd(double x, double D, double d2) const {
		return 2.0 * D / (_Re * U(x) / _U_inf * d2);
	}

	double BoundaryLayer::H12(double h32) const {

		if (h32 >= 1.57258) {
			return 79.870845 - 89.582142 * h32 + 25.715786 * h32 * h32;
		} else if (h32 >= 1.51509) {
			return 4.02922 - (583.60182 - 724.55916 * h32
					  + 227.18220 * h32 * h32)
					  * std::sqrt(h32 - 1.51509);
		} else {
			return 0.0;	// !!!!????
		}

	}

	double BoundaryLayer::Eps(double h32) const {

		double h12 = H12(h32);

		if (H32 >= 1.57258) {
			return 1.372391 - 4.226253 * h32 + 2.221687 * h32 * h32;
		} else if (h32 >= 1.51509) {
			return 2.512589 - 1.686095 * h12 + 0.391541 * h12 * h12
					  - 0.031729 * h12 * h12 * h12;
		} else {
			return 0.0;	// !!!!????
		}
		
	}

	double BoundaryLater::D(double H32) const {
		return 7.853976 - 10.260551 * H32 + 3.418898 * H32 * H32;
	}

}


