/**
 * File:   Airfoil.hpp
 * Author: kenobi
 *
 * Created on July 23, 2014, 12:51 PM
 */

#ifndef INCLUDE_AIRFOIL_HPP
#define INCLUDE_AIRFOIL_HPP

#include "Interpolator.hpp"
#include "NacaAirfoil.hpp"
#include "Panel.hpp"
#include "ThinAirfoil.hpp"
#include "Vector.hpp"
#include <cmath>
#include <string>

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
	class Airfoil {

		enum CoordsType {
			naca = 1,
			points = 2
		};

	public:
		Airfoil(const NacaAirfoil& naca_airfoil);
		Airfoil(const std::string& data_airfoil_file, double chord);
		~Airfoil();

		Vector<Panel2D > getPanels(double density) const;
		double getChord() const;
		double getPointBottom(double x) const;
		double getPointTop(double x) const;

		void writeCoordsToFile(const std::string& file_name,
									unsigned int n) const;

	private:
		CoordsType _coords_type;
		NacaAirfoil *_naca_airfoil;
		Interpolator1D *_interp_top;
		Interpolator1D *_interp_bottom;

		double _chord;
		
	};

}

#endif	/* INCLUDE_AIRFOIL_HPP */

