/**
 * File:   Airfoil.cpp
 * Author: kenobi
 *
 * Created on July 23, 2014, 12:51 PM
 */

#include "Airfoil.hpp"
#include <cmath>
#include <fstream>
#include <sstream>

namespace nde {

	Airfoil::Airfoil(const NacaAirfoil& naca_airfoil) {
		_coords_type = Airfoil::naca;
		_naca_airfoil = new NacaAirfoil(naca_airfoil);
		_chord = _naca_airfoil->getChord();
	}

	Airfoil::Airfoil(const std::string& data_airfoil_file, double chord) :
						_chord(chord) {

		_coords_type = Airfoil::points;

		// read input points
		std::ifstream in_file(data_airfoil_file.c_str(), std::ifstream::in);
		unsigned int n_in = 0;
		Vector<double> x_all(1000), y_all(1000);
		while (in_file.good()) {
			std::string line;	std::getline(in_file, line);
			std::istringstream iss(line);
			iss >> x_all(n_in) >> y_all(n_in);
			n_in += 1;
		}
		n_in -= 1;
		in_file.close();

		// initialize interpolating splines
		unsigned int n = n_in / 2 + 1;
		Vector<double> x_bottom(n), y_bottom(n);
		Vector<double> x_top(n), y_top(n);

		unsigned int k = n_in / 2;

		for (unsigned int i = 0; i < n; ++i) {
			x_bottom(i) = x_all(n - 1 - i);
			y_bottom(i) = y_all(n - 1 - i);
			x_top(i) = x_all(k + i);
			y_top(i) = y_all(k + i);
		}

		_interp_top = new Interpolator1D(x_top, y_top, SPLINE_MONOTONE);
		_interp_bottom = new Interpolator1D(x_bottom, y_bottom, SPLINE_MONOTONE);

	}

	Airfoil::~Airfoil() {
		
		if (_coords_type == Airfoil::naca) {
			delete _naca_airfoil;
		} else if (_coords_type == Airfoil::points) {
			delete _interp_top;
			delete _interp_bottom;
		}

	}


	Vector<Panel2D> Airfoil::getPanels(double density) const {

		int num_panels = floor(1.0 / density);
		double dx = _chord / double(num_panels);

		Vector<Panel2D > x(2 * num_panels);

		for (int i = 0; i < num_panels; ++i) {

			Vector<double> p_start_bottom(2);
			Vector<double> p_end_bottom(2);
			p_start_bottom(0) = density * (num_panels - i) * _chord;
			p_start_bottom(1) = getPointBottom(p_start_bottom(0));
			p_end_bottom(0) = density * (num_panels - i - 1) * _chord;
			p_end_bottom(1) = getPointBottom(p_end_bottom(0));

			x(i).setPoints(p_start_bottom, p_end_bottom);

			Vector<double> p_start_top(2);
			Vector<double> p_end_top(2);
			p_start_top(0) = density * i * _chord;
			p_start_top(1) = getPointTop(p_start_top(0));
			p_end_top(0) = density * (i + 1) * _chord;
			p_end_top(1) = getPointTop(p_end_top(0));

			x(num_panels+i).setPoints(p_start_top, p_end_top);

		}

		return x;

	}

	double Airfoil::getChord() const {
		return _chord;
	}

	double Airfoil::getPointBottom(double x) const {
		if (_coords_type == Airfoil::naca) {
			return _naca_airfoil->bottom(x);
		} else if (_coords_type == Airfoil::points) {
			return (*_interp_bottom)(x / _chord) * _chord;
		}
	}

	double Airfoil::getPointTop(double x) const {
		if (_coords_type == Airfoil::naca) {
			return _naca_airfoil->top(x);
		} else if (_coords_type == Airfoil::points) {
			return (*_interp_top)(x / _chord) * _chord;
		}
	}

	void Airfoil::writeCoordsToFile(const std::string& file_name,
											 unsigned int n) const {

		std::ofstream out_file(file_name.c_str(), std::ofstream::out);

		// bottom surface
		for (unsigned int i = 0; i < n; ++i) {
			double x = _chord * (1.0 - double(i)/double(n-1));
			out_file << x << " " << getPointBottom(x) << std::endl;
		}

		// top surface
		for (unsigned int i = 0; i < n; ++i) {
			double x = _chord * double(i)/double(n-1);
			out_file << x << " " << getPointTop(x) << std::endl;
		}

		out_file.close();
		
	}

} /* end of namespace mrm */

