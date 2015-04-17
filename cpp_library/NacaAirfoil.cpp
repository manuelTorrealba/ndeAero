/**
 * File:   NacaAirfoil.cpp
 * Author: kenobi
 *
 * Created on Apr 17, 2015, 11:04 AM
 */

#include "NacaAirfoil.hpp"
#include <cmath>

namespace nde {


	NacaAirfoil::NacaAirfoil(double chord, const Vector<int>& naca_codenum)
	  : _chord(chord_in), _naca_codenum(naca_codenum_in) {

		_num_digits = naca_codenum.size();

		//
		// calculate maximum thickness
		//
		if (_num_digits == 4) // for the four digits naca airfoil
									// digits 3 and 4 represent maximum thickness
			_max_thickness = (double(naca_codenum(2)) * 10.
									+ double(naca_codenum(3))) / 100.;
		else if (_num_digits == 5) // for the five digits naca airfoil
											// digits 4 and 5 represent maximum thickness
			_max_thickness = (double(naca_codenum(3)) * 10.
									+ double(naca_codenum(4))) / 100.;

		//
		// calculate the maximum camber value and position
		//
		if (_num_digits == 4) {
		// four digits naca airfoil

			// maximum value of camber line
			_max_camber_val = double(_naca_codenum(0)) / 100.;
			// maximum camber position
			_max_camber_x = double(_naca_codenum(1)) / 10.;

		} else if (_num_digits == 5) {
		//five digits naca airfoil

			int q = int(_naca_codenum(2));

			if (q==0) {
			//standard camber line

				helperCtor5digits(_max_camber_val, _max_camber_x);

			} else if (q==1) {
			//reflexed camber line -> achieves small pitch

				helperCtor5digitsReflex(_max_camber_val, _max_camber_x);

			}

		}

	}

	void helperCtor5digits(double& max_camber_val, double& max_camber_x) const {

		double l = double(_naca_codenum(0)) / 10.;
		double p = double(_naca_codenum(1)) / 10.;
		double cl_ideal = 1.5 * l; // design lift coefficient
		double xf = 0.5 * p; // maximum camber position

		// maximum camber position designator calculation
		// by fixed point iterations
		double m = xf;
		double err;
		do {
			double m_new = xf / (1.0 - std::sqrt(m/3.0));
			err = std::abs(m - m_new);
			m = m_new;
		} while (err > 1.e-6);

		max_camber_x = m;

		// _max_camber_val quantity chosen to avoid leading edge singularity
		double q1 = (3.0 * m - 7.0 * m2  + 8.0 * m3 - 4.0 * m4)
					 / std::sqrt(m - m2);
		double q2 = -1.5 * ( 1.0 - 2.0 * m) *
						(0.5 * M_PI - std::asin(1.0 - 2.0 * m));
		max_camber_val = cl_ideal / (q1 + q2);

	}

	void helperCtor5digitsReflex(double& max_camber_val,
									double& max_camber_x) const {

		double l = double(_naca_codenum(0)) / 10.;
		double p = double(_naca_codenum(1)) / 10.;
		double cl_ideal = 1.5 * l; // design lift coefficient
		double xf = 0.5 * p; // maximum camber position

		double m; // m is found to give Cm(c/4) = 0 from thin airfoil theory.
		double k1; // k1 is calculated to give Cl(ideal) = design lift coeffcient

	}

	double NacaAirfoil::camber(double x) const {

		if (_num_digits == 4)
		//four digits naca airfoil

			return camber4digits(x);

		else if (_num_digits == 5) {
		//five digits naca airfoil

			int q = int(_naca_codenum(2));

			if (q==0) {
			//standard camber line

				return camber5digits(x);

			} else if (q==1) {
			//reflexed camber line -> achieves small pitch

				return camber5digitsReflex(x);
			}

		} else

			return 0.0;

	}

	double NacaAirfoil::camber4digits(double x) const {

		//four digits naca airfoil camber -> easy made!

		// short names
		double m = _max_camber_val;
		double p = _max_camber_x;

		double y = 0.0;
		if (x < p * _chord)
			y = m * x / (p * p)*(2.0 * p - x / _chord);
		else
			y = m * (_chord - x) / ((1.0 - p) * (1.0 - p))
			    * (1.0 + x / _chord - 2.0 * p);

		return y;

	}

	double NacaAirfoil::dcamberdx4digits(double x) const {


	}

	double NacaAirfoil::camber5digits(double x) const {

		// short names
		double m = _max_camber_x;
		double k1 = _max_camber_val;

		// different powers of _max_camber_x needed
		double m2 = m * m;
		double m3 = m2 * m;
		double m4 = m3 * m;


		// all camber parameters defined -> calculate camber!
		double t = x / _chord;
		double t2 = t * t;
		double t3 = t2 * t;

		double y;
		if (t < m) {
			y = k1 * (t3 - 3.0 * m * t2 + m2 * (3.0 - m) * t);
		} else {
			y = k1 * m3 * (1.0 - t);
		}

		return _chord * y;

	}

	double NacaAirfoil::camber5digitsReflex(double x) const {

		// short names
		double m = _max_camber_x;
		double k1 = _max_camber_val;

		double p = double(_naca_codenum(1)) / 10.;
		double xf = 0.5 * p; // maximum camber position
		double k2 = 3.0 * std::pow(m - xf, 2);

		// camber line equation
		double y = 0.0;
		double t = x / _chord;
		if (t <= m)
			y = std::pow(t - m, 3) - k2 * t + std::pow(m, 3);
		else
			y = k2 * std::pow(t - m, 3) - k2 * t + std::pow(m, 3);

		return y * _chord * k1 / 6.0;

	}

	double NacaAirfoil::dcamberdx5digits(double x) const;
	double NacaAirfoil::dcamberdx5digitsReflex(double x) const;

	double NacaAirfoil::thickness(double x) const {

		double t = x / chord;
		double y = _max_thickness / 0.2 * (0.2969 * std::sqrt(t) - 0.1260 * t
		        - 0.3516 * t * t + 0.2843 * t * t * t
		        - 0.1036 * t * t * t * t);
		//last coefficient modified from -0.1015 to -0.1036
		//to get zero thickness at the trailing edge while
		//modifying the shape of the airfoil as little
		//as possible.
		return y * chord;
	}

} /* end of namespace nde */

