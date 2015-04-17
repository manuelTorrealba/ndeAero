#include "ThinAirfoil.hpp"
#include <cmath>

namespace nde {

	ThinAirfoil::ThinAirfoil(unsigned int n) {
		_n_steps = std::pow(2, n);
	}

	double ThinAirfoil::calcCL(double angle_attack) const {
		return 2.0 * M_PI * (angle_attack + calcA(0) + 0.5 * calcA(1));
	}

	double ThinAirfoil::calcCM(double angle_attack, double tm) const {
		double cm0 = - 0.5 * M_PI * (angle_attack + calcA(0)
											+ calcA(1) - 0.5 * calcA(2));
		double cl = calcCL(angle_attack);
		return cm0 + cl * tm;
	}

	double ThinAirfoil::calcA(unsigned int n) const {
		
		double dt = M_PI / _n_steps;
	
		double a = 0.0;
		for (unsigned int i = 0; i < _n_steps; ++i) {
			
			double t = dt*double(i);
			double c = std::cos(t*double(n));
			double d = dCamberDx(t);
			a += c * d;
		}

		a /= M_PI;

		if (n == 0)
			a *= -1.0;
		else
			a *= 2.0;

		return a;

	}

}

