/**
 * File:   NacaAirfoil.hpp
 * Author: kenobi
 *
 * Created on Apr 17, 2015, 11:02 AM
 */

#ifndef INCLUDE_NACA_AIRFOIL_HPP
#define INCLUDE_NACA_AIRFOIL_HPP

#include "Vector.hpp"

namespace nde {

	/**
	  * Naca Airfoil structure
	  */
	struct NacaAirfoil {
	public:
      NacaAirfoil(double chord, const Vector<int>& naca_codenum);

		double getChord() const {
			return _chord;
		}

      double top(double x) const {
      	return camber(x) + thickness(x);
      }

      double bottom(double x) const {
      	return camber(x) - thickness(x);
		}

	protected:
		double camber(double x) const;
		double thickness(double x) const;
		double dcamberdx(double x) const;

	private:
   	double _chord;
		Vector<int> _naca_codenum;

		unsigned int _num_digits;
		double _max_thickness;
		double _max_camber_val;
		double _max_camber_x; // maximum camber position

		// TODO:
		// two more structures
		// Naca4Digits
		// Naca5Digits
		double camber4digits(double x) const;
		double camber5digits(double x) const;
		double camber5digitsReflex(double x) const;
		double dcamberdx4digits(double x) const;
		double dcamberdx5digits(double x) const;
		double dcamberdx5digitsReflex(double x) const;

		void helperCtor5digits(double& max_camber_val,
									double & max_camber_x) const;
		void helperCtor5digitsReflex(double& max_camber_val,
										double& max_camber_x) const;

	};

} /* end of namespace nde */

