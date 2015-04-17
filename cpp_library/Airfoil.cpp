/**
 * File:   Airfoil.cpp
 * Author: kenobi
 *
 * Created on July 23, 2014, 12:51 PM
 */

#include "Airfoil.hpp"
#include <cmath>

namespace nde {

	Airfoil::Airfoil(const NacaAirfoil& naca_airfoil_in)
						: ThinAirfoil(5), naca_airfoil(naca_airfoil_in) {
	  ;
	}

	Vector<Panel2D > Airfoil::getPanels(double density) const {

		int num_panels = floor(1.0 / density);
		double dx = naca_airfoil.chord / double(num_panels);

		Vector<Panel2D > x(2 * num_panels);

		for (int i = 0; i < num_panels; ++i) {
			Vector<double> p_start_bottom(2);
			Vector<double> p_end_bottom(2);
			p_start_bottom(0) = density * (num_panels - i) * naca_airfoil.chord;
			p_start_bottom(1) = naca_airfoil.bottom(p_start_bottom(0));
			p_end_bottom(0) = density * (num_panels - i - 1) * naca_airfoil.chord;
			p_end_bottom(1) = naca_airfoil.bottom(p_end_bottom(0));

			x(i).setPoints(p_start_bottom, p_end_bottom);

			Vector<double> p_start_top(2);
			Vector<double> p_end_top(2);
			p_start_top(0) = density * i * naca_airfoil.chord;
			p_start_top(1) = naca_airfoil.top(p_start_top(0));
			p_end_top(0) = density * (i + 1) * naca_airfoil.chord;
			p_end_top(1) = naca_airfoil.top(p_end_top(0));

			x(num_panels+i).setPoints(p_start_top, p_end_top);
		}

		return x;

	}

}
