/**
 * File:   ThinAirfoil.hpp
 * Author: kenobi
 *
 * Created on April 16, 2015, 12:51 PM
 */

#ifndef INCLUDE_THIN_AIRFOIL_HPP
#define INCLUDE_THIN_AIRFOIL_HPP

namespace nde {

	/**
	  * This class calculates results from Thin Airfoil Theory.
	  * Integrals are calculated numerically.
	  */
	class ThinAirfoil {
	public:
		ThinAirfoil(unsigned int n);
		double calcCL(double angle_attack) const;
		double calcCM(double angle_attack, double tm) const;

		virtual double camber(double t) const = 0;
		virtual	double dCamberDx(double t) const = 0;


	protected:
		unsigned int _n_steps;
		double calcA(unsigned int n) const;
	};

} /* end of nde */

#endif

