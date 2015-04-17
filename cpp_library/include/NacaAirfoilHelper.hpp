/**
 * File:   NacaAirfoilHelper.hpp
 * Author: kenobi
 *
 * Created on April 16, 2015, 12:51 PM
 */

#ifndef INCLUDE_NACA_AIRFOIL_HELPER_HPP
#define INCLUDE_NACA_AIRFOIL_HELPER_HPP

#include"Vector.hpp"
#include"Panel.hpp"
#include<cmath>

namespace nde {

	/**
	  * Naca
	  */
	class NacaFiveDigitsOptCamber {
	public:
		NacaFiveDigitsOptCamber(double xf);
		double calcThinAirfoilMoment
	private:
		double _xf;

	}

}

#endif

