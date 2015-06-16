#include "AerodynamicBody2D.hpp"
#include "Airfoil.hpp"
#include "AirfoilBoundaryLayer.hpp"
#include "Test.hpp"

#include <iostream>

using namespace std;

namespace nde {

bool testNACAAirfoil(const Vector<unsigned int> &naca_codenum,
							const std::string& out_file_name) {

	cout << "/******************************************************/" << endl;
	cout << "/* NACA airfoil test" << endl;
	cout << "/******************************************************/" << endl;

	cout << "Build the Airfoil and output coordinates in "
		  << out_file_name << "_coords.dat" << endl;

	std::string out_airfoil(out_file_name + "_coords.dat");

	nde::NacaAirfoil naca_airfoil(1.0, naca_codenum);
	nde::Airfoil airfoil(naca_airfoil);
	airfoil.writeCoordsToFile(out_airfoil,50);


	// Panelize the airfoil
	nde::Vector<nde::Panel2D> panels = airfoil.getPanels(0.1);
	nde::AerodynamicBody2D body2D(1.0, panels);

	cout << "AoA, Cl1, Cl2, Cl3, M1, M2, M3, delta_Cd, delta_Cl" << endl;

	// test the solution for angles of attack between -10 and +10 degrees
	for (int a = -20; a<= 20; ++a) {

		double a_rad = double(a) * M_PI / 180.0;
		body2D.calcPotentialFlow(a_rad,
										 nde::DIRICHLET_CONSTANT_DOUBLETS);
		nde::Vector<double> F1 = body2D.getForceCoeffs();
		double M1 = body2D.getMomentCoeff();

		body2D.calcPotentialFlow(a_rad,
										 nde::DIRICHLET_CONSTANT_SOURCES_AND_DOUBLETS);
		nde::Vector<double> F2 = body2D.getForceCoeffs();
		double M2 = body2D.getMomentCoeff();

		body2D.calcPotentialFlow(a_rad,
										 nde::NEUMANN_CONSTANT_SOURCES_AND_VORTEX);
		nde::Vector<double> F3 = body2D.getForceCoeffs();
		double M3 = body2D.getMomentCoeff();

		/*********************************************************************/
		/* boundary layer
		/*********************************************************************/
		nde::AirfoilBoundaryLayer airfoil_boundary_layer(100. / 3.6, // 100 km/h
														2.0, // chord = 2 m
														body2D);

		airfoil_boundary_layer.solveBoundaryLayerFlow();
		double delta_Cd = airfoil_boundary_layer.getDeltaCD();
		double delta_Cl = airfoil_boundary_layer.getDeltaCL();

		cout << a << ", " << F1(1) << ", " << F2(1) << ", " << F3(1)
			  << ", " << M1 << ", " << M2 << ", " << M3 << ", "
			  << delta_Cd << ", " << delta_Cl << endl;

	}

	/*
		cout << "Cl Dirichlet 1 = (" << F1(0) << ", " << F1(1) << ")" << endl;
		cout << "Cm Dirichlet 1 = (" << M1 << ")" << endl;

		body2D.writeResultsToFile("NACA0012_results_dd.dat");

		cout << "Cl Dirichlet 2 = (" << F2(0) << "," << F2(1) << ")" << endl;
		cout << "Cm Dirichlet 2 = (" << M2 << ")" << endl;

		body2D.writeResultsToFile("NACA0012_results_dsd.dat");

		cout << "Cl Neumann = (" << F3(0) << "," << F3(1) << ")" << endl;
		cout << "Cm Neumann = (" << M3 << ")" << endl;

		body2D.writeResultsToFile("NACA0012_results_nsv.dat");
	*/

	return true;

}

} /* end of namespace nde */

