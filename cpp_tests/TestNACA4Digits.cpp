#include "AerodynamicBody2D.hpp"
#include "Airfoil.hpp"
#include "BoundaryLayer.hpp"
#include "Test.hpp"

#include <iostream>

using namespace std;

namespace nde {

bool testNACA4Digits() {

	cout << "/******************************************************/" << endl;
	cout << "/* NACA 4-digits airfoil test" << endl;
	cout << "/******************************************************/" << endl;


	cout << "Build the NACA0012 Airfoil \
				-> output coordinates in NACA0012.dat" << endl;

	nde::Vector<unsigned int> naca_codenum(4);
	naca_codenum(0) = 0;	naca_codenum(1) = 0;
	naca_codenum(2) = 1; naca_codenum(3) = 2;
	nde::NacaAirfoil naca_airfoil(1.0, naca_codenum);
	nde::Airfoil airfoil(naca_airfoil);
	airfoil.writeCoordsToFile("NACA0012.dat",50);

	cout << "Panelize. 3 degrees angle of attack" << endl;
	nde::Vector<nde::Panel2D> panels = airfoil.getPanels(0.01);
	nde::AerodynamicBody2D body2D(1.0, panels);

	body2D.calcPotentialFlow(3.0 * M_PI / 180.0,
									 nde::DIRICHLET_CONSTANT_DOUBLETS);
	nde::Vector<double> F1 = body2D.getForceCoeffs();
	cout << "Force Dirichlet 1 = (" << F1(0) << ", " << F1(1) << ")" << endl;

	body2D.writeResultsToFile("NACA0012_results_dd.dat");

	body2D.calcPotentialFlow(3.0 * M_PI / 180.0,
									 nde::DIRICHLET_CONSTANT_SOURCES_AND_DOUBLETS);
	nde::Vector<double> F2 = body2D.getForceCoeffs();
	cout << "Force Dirichlet 2 = (" << F2(0) << "," << F2(1) << ")" << endl;

	body2D.writeResultsToFile("NACA0012_results_dsd.dat");

	body2D.calcPotentialFlow(3.0 * M_PI / 180.0,
									 nde::NEUMANN_CONSTANT_SOURCES_AND_VORTEX);
	nde::Vector<double> F3 = body2D.getForceCoeffs();
	cout << "Force Neumann = (" << F3(0) << "," << F3(1) << ")" << endl;

	body2D.writeResultsToFile("NACA0012_results_nsv.dat");

	/*********************************************************************/
	/* boundary layer
	/*********************************************************************/
	nde::AirfoilBoundaryLayer airfoil_boundary_layer(100. / 3.6, // 100 km/h
													2.0, // chord = 2 m
													body2D.getControlPointsDistances(),
													body2D.getVelocities());


	return true;

}

} /* end of namespace nde */

