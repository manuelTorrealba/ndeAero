/**
 * File:   Airfoil.hpp
 * Author: kenobi
 *
 * Created on July 23, 2014, 12:51 PM
 */

#ifndef INCLUDE_AIRFOIL_HPP
#define INCLUDE_AIRFOIL_HPP

#include "Panel.hpp"
#include "NacaAirfoil.hpp"
#include "ThinAirfoil.hpp"
#include "Vector.hpp"
#include <cmath>

namespace nde {

	/**
	  * Generic Airfoil class
	  * This routine panelizes a given airfoil curve.
	  * It can be initialized with a Naca Airfoil.
	  *
	  * TODO
	  *	- initialize with a generic coordinates vector.
	  *	- define flap surfaces and flap deflexion.
     *
	  */
	class Airfoil : public ThinAirfoil {
	public:
		Airfoil(const NacaAirfoil& naca_airfoil_in);
		Vector<Panel2D > getPanels(double density) const;
		double getChord() const;
		double getPointBottom(double x) const;
		double getPointTop(double x) const;

	protected:
		// inherited from ThinAirfoil
		virtual	double dCamberDx(double t) const;

	private:
		NacaAirfoil naca_airfoil;
	};

}

#endif	/* INCLUDE_AIRFOIL_HPP */

