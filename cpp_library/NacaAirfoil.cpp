/**
 * File:   NacaAirfoil.cpp
 * Author: kenobi
 *
 * Created on Apr 17, 2015, 11:04 AM
 */

#include "EquationSolver.hpp"
#include "NacaAirfoil.hpp"
#include <cmath>

namespace nde {


	NacaAirfoil::NacaAirfoil(double chord,
								const Vector<unsigned int>& naca_codenum)
	  : _chord(chord), _naca_codenum(naca_codenum) {

		helperCtor();

	}


	NacaAirfoil::NacaAirfoil(const NacaAirfoil& naca_airfoil) {
		
		this->_chord = naca_airfoil._chord;
		this->_naca_codenum = naca_airfoil._naca_codenum;
		helperCtor();

	}

	void NacaAirfoil::helperCtor() {

		_num_digits = _naca_codenum.size();

		if (_num_digits == 4) {

			_max_thickness = double(_naca_codenum(2) * 10 + _naca_codenum(3)) / 100.;
			_naca_thin_airfoil = new NacaFourDigits(_naca_codenum);

		} else if (_num_digits == 5) {

			_max_thickness = double(_naca_codenum(3) * 10 + _naca_codenum(4)) / 100.;

			int q = int(_naca_codenum(2));
			if (q == 0) // standard camber line
				_naca_thin_airfoil = new NacaFiveDigits(_naca_codenum);
			else if (q == 1) // reflexed camber line
				_naca_thin_airfoil = new NacaFiveDigitsReflx(_naca_codenum);

		}

	}

	NacaAirfoil::~NacaAirfoil() {

		delete _naca_thin_airfoil;

	}


	double NacaAirfoil::camber(double x) const {

			return _chord * _naca_thin_airfoil->camber(x / _chord);

	}


	double NacaAirfoil::thickness(double x) const {

		double t = x / _chord;
		double y = _max_thickness / 0.2 * (0.2969 * std::sqrt(t) - 0.1260 * t
		        - 0.3516 * t * t + 0.2843 * t * t * t
		        - 0.1036 * t * t * t * t);
		//last coefficient modified from -0.1015 to -0.1036
		//to get zero thickness at the trailing edge while
		//modifying the shape of the airfoil as little
		//as possible.
		return y * _chord;
	}


	/**
	  * NACA 4-digits airfoil
	  */
	NacaFourDigits::NacaFourDigits(Vector<unsigned int> digits)
											:ThinAirfoil(14), _digits(digits) {

		//
		// calculate the maximum camber value and position
		//
		_max_camber_val = double(_digits(0)) / 100.;
		_max_camber_x = double(_digits(1)) / 10.;

	}


	double NacaFourDigits::camber(double t) const {
		
		//four digits naca airfoil camber -> easy made!

		// short names
		double m = _max_camber_val;
		double p = _max_camber_x;

		double y = 0.0;
		if (t < p)
			y = m * t / (p * p) * (2.0 * p - t);
		else
			y = m * (1.0 - t) / ((1.0 - p) * (1.0 - p)) * (1.0 + t - 2.0 * p);

		return y;

	}

	double NacaFourDigits::dCamberDx(double t) const {

		// short names
		double m = _max_camber_val;
		double p = _max_camber_x;

		double y = 0.0;
		if (t < p)
			y = 2.0 * m / (p * p) * (p - t);
		else
			y = 2.0 * m * ((1.0 - p) * (1.0 - p)) * (p - t);

		return y;

	}


	
	/**
	  * NACA 5-digits airfoil non-reflexed
	  */
	NacaFiveDigits::NacaFiveDigits(Vector<unsigned int> digits)
											: ThinAirfoil(14), _digits(digits) {

		helperCtor(_max_camber_val, _max_camber_x);

		std::cout << _max_camber_val << "," << _max_camber_x << std::endl;

	}

	void NacaFiveDigits::helperCtor(double& max_camber_val,
											double& max_camber_x) const {

		double l = double(_digits(0)) / 10.;
		double p = double(_digits(1)) / 10.;
		double cl_ideal = 1.5 * l; // design lift coefficient
		double xf = 0.5 * p; // maximum camber position designator

		// maximum camber position calculation
		// by fixed point iterations
		double m = xf;
		double err;
		do {
			double m_new = xf / (1.0 - std::sqrt(m/3.0));
			err = std::abs(m - m_new);
			m = m_new;
		} while (err > 1.e-6);

		max_camber_x = m;

		double m2 = m * m;
		double m3 = m2 * m;
		double m4 = m3 * m;

		// _max_camber_val quantity chosen to avoid leading edge singularity
		double q1 = (3.0 * m - 7.0 * m2  + 8.0 * m3 - 4.0 * m4)
					 / std::sqrt(m - m2);
		double q2 = -1.5 * ( 1.0 - 2.0 * m) *
						(0.5 * M_PI - std::asin(1.0 - 2.0 * m));
		max_camber_val = cl_ideal / (q1 + q2);

	}

	double NacaFiveDigits::camber(double t) const {

		// short names
		double m = _max_camber_x;
		double k1 = _max_camber_val;

		// different powers of _max_camber_x needed
		double m2 = m * m;
		double m3 = m2 * m;
		double m4 = m3 * m;


		// all camber parameters defined -> calculate camber!
		double t2 = t * t;
		double t3 = t2 * t;

		double y;
		if (t < m) {
			y = k1 * (t3 - 3.0 * m * t2 + m2 * (3.0 - m) * t);
		} else {
			y = k1 * m3 * (1.0 - t);
		}

		return y;

	}


	double NacaFiveDigits::dCamberDx(double t) const {

		// short names
		double m = _max_camber_x;
		double k1 = _max_camber_val;

		// different powers of _max_camber_x needed
		double m2 = m * m;
		double m3 = m2 * m;
		double m4 = m3 * m;


		// all camber parameters defined -> calculate camber!
		double t2 = t * t;

		double y;
		if (t < m) {
			y = k1 * (3.0 * t2 - 6.0 * m * t + m2 * (3.0 - m));
		} else {
			y = -k1 * m3;
		}

		return y;

	}



	/**
	  * NACA 5-digits airfoil reflexed
	  */
	NacaFiveDigitsReflx::NacaFiveDigitsReflx(Vector<unsigned int> digits)
														: ThinAirfoil(14), _digits(digits) {

		// _max_camber_x (-> max camber position) is found to give Cm(c/4) = 0
		// from thin airfoil theory.
		double p = double(_digits(1)) / 10.;
		double xf = 0.5 * p; // maximum camber position designator
									// search for the root which is on the right of xf!

		_solver_objective = 1;
		_max_camber_x = secant(std::min(1.2 * xf, 0.9),
									  std::min(2.0 * xf, 1.0));

		// _max_camber_val is calculated to give
		// Cl(ideal) = design lift coeffcient
		_solver_objective = 2;
		_max_camber_val = secant(10, 20);

		std::cout << _max_camber_val << "," << _max_camber_x << std::endl;

	}

	double NacaFiveDigitsReflx::solverF(double x) {

		if (_solver_objective == 1) {
			// _max_camber_val is irrelevant for equation Cm(c/4) = 0
			_max_camber_val = 1.0;
			_max_camber_x = x;
			return calcCM(0.0, 0.25); // angle of attack parameter is irrelevant for
												// this problem
		} else if (_solver_objective == 2) {

			// _max_camber_val is calculated to give
			// Cl(ideal) = design lift coeffcient
			double l = double(_digits(0)) / 10.;
			double cl_ideal = 1.5 * l; // design lift coefficient
			_max_camber_val = x;
			return M_PI * calcA(1) - cl_ideal;

		} else {

			//TODO:			std::runtimeerror("");

		}

	}

	double NacaFiveDigitsReflx::camber(double t) const {

		// short names
		double m = _max_camber_x;
		double k1 = _max_camber_val;
		double m3 = std::pow(m, 3);

		double p = double(_digits(1)) / 10.;
		double xf = 0.5 * p; // maximum camber position designator
		double k2 = (3.0 * std::pow(m - xf, 2) - m3) / std::pow(1.0 - m, 3);

		// camber line equation
		double y = 0.0;
		if (t <= m)
			y = std::pow(t - m, 3) - k2 * std::pow(1.0 - m, 3) * t - m3 * t + m3;
		else
			y = k2 * std::pow(t - m, 3) - k2 * std::pow(1.0 - m, 3) * t
			  - m3 * t + m3;

		return y * k1 / 6.0;

	}

	double NacaFiveDigitsReflx::dCamberDx(double t) const {

		// short names
		double m = _max_camber_x;
		double k1 = _max_camber_val;
		double m3 = std::pow(m, 3);

		double p = double(_digits(1)) / 10.;
		double xf = 0.5 * p; // maximum camber position designator
		double k2 = (3.0 * std::pow(m - xf, 2) - m3) / std::pow(1.0 - m, 3);

		// camber line equation
		double y = 0.0;
		if (t <= m)
			y = 3.0 * std::pow(t - m, 2) - k2 * std::pow(1.0 - m, 3) - m3;
		else
			y = k2 * 3.0 * std::pow(t - m, 2) - k2 * std::pow(1.0 - m, 3) - m3;

		return y * k1 / 6.0;

	}

} /* end of namespace nde */

