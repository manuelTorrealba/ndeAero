#include "Test.hpp"
#include <iostream>

/* nde: numerical design engineering */

using namespace std;

int main(int narg, char** arg) {


	// test thin airfoil theory formulas
	cout << nde::testThinAirfoilTheory() << endl;

	// test the ODE solver
	cout << nde::testODESolver() << endl;

	// test the NACA 4 digits series
	cout << nde::testNACA4Digits() << endl;


//	{

//	cout << "/******************************************************/" << endl;
//	cout << "/* NACA 5-digits airfoil test" << endl;
//	cout << "/******************************************************/" << endl;

//	nde::Vector<unsigned int> naca_codenum(5);
//	naca_codenum(0) = 2;	naca_codenum(1) = 3;
//	naca_codenum(2) = 0; naca_codenum(3) = 1; naca_codenum(4) = 2;
//	nde::NacaAirfoil naca_airfoil(1.0, naca_codenum);
//	nde::Airfoil airfoil(naca_airfoil);
//	airfoil.writeCoordsToFile("NACA23012.dat",50);

//	nde::Vector<nde::Panel2D> panels = airfoil.getPanels(0.01);
//	nde::AerodynamicBody2D body2D(1.0, panels);


//	cout << "Angle of attack = 0.0" << endl;

//	body2D.calcPotentialFlow(0.0, nde::DIRICHLET_CONSTANT_DOUBLETS);
//	nde::Vector<double> F1 = body2D.getForceCoeffs();
//	cout << "Force Dirichlet 1 = (" << F1(0) << ", " << F1(1) << ")" << endl;

//	body2D.calcPotentialFlow(0.0, nde::DIRICHLET_CONSTANT_SOURCES_AND_DOUBLETS);
//	nde::Vector<double> F2 = body2D.getForceCoeffs();
//	cout << "Force Dirichlet 2 = (" << F2(0) << "," << F2(1) << ")" << endl;

//	body2D.calcPotentialFlow(0.0, nde::NEUMANN_CONSTANT_SOURCES_AND_VORTEX);
//	nde::Vector<double> F3 = body2D.getForceCoeffs();
//	cout << "Force Neumann = (" << F3(0) << "," << F3(1) << ")" << endl;


//	cout << "Angle of attack = 0.0873" << endl;

//	body2D.calcPotentialFlow(5.0 * M_PI / 180.0,
//									 nde::DIRICHLET_CONSTANT_DOUBLETS);
//	nde::Vector<double> F1_1 = body2D.getForceCoeffs();
//	cout << "Force Dirichlet 1 = (" << F1_1(0) << ", " << F1_1(1) << ")" << endl;

//	body2D.calcPotentialFlow(5.0 * M_PI / 180.0,
//									 nde::DIRICHLET_CONSTANT_SOURCES_AND_DOUBLETS);
//	nde::Vector<double> F2_1 = body2D.getForceCoeffs();
//	cout << "Force Dirichlet 2 = (" << F2_1(0) << "," << F2_1(1) << ")" << endl;

//	body2D.calcPotentialFlow(5.0 * M_PI / 180.0,
//									 nde::NEUMANN_CONSTANT_SOURCES_AND_VORTEX);
//	nde::Vector<double> F3_1 = body2D.getForceCoeffs();
//	cout << "Force Neumann = (" << F3_1(0) << "," << F3_1(1) << ")" << endl;

//	}


//	{

//	cout << "/******************************************************/" << endl;
//	cout << "/* NACA 5-digits reflexed airfoil test" << endl;
//	cout << "/******************************************************/" << endl;

//	nde::Vector<unsigned int> naca_codenum(5);
//	naca_codenum(0) = 2;	naca_codenum(1) = 3;
//	naca_codenum(2) = 1; naca_codenum(3) = 1; naca_codenum(4) = 2;
//	nde::NacaAirfoil naca_airfoil(1.0, naca_codenum);
//	nde::Airfoil airfoil(naca_airfoil);
//	airfoil.writeCoordsToFile("NACA23112.dat",50);

//	nde::Vector<nde::Panel2D> panels = airfoil.getPanels(0.01);
//	nde::AerodynamicBody2D body2D(1.0, panels);


//	cout << "Angle of attack = 0.0" << endl;

//	body2D.calcPotentialFlow(0.0, nde::DIRICHLET_CONSTANT_DOUBLETS);
//	nde::Vector<double> F1 = body2D.getForceCoeffs();
//	cout << "Force Dirichlet 1 = (" << F1(0) << ", " << F1(1) << ")" << endl;

//	body2D.calcPotentialFlow(0.0, nde::DIRICHLET_CONSTANT_SOURCES_AND_DOUBLETS);
//	nde::Vector<double> F2 = body2D.getForceCoeffs();
//	cout << "Force Dirichlet 2 = (" << F2(0) << "," << F2(1) << ")" << endl;

//	body2D.calcPotentialFlow(0.0, nde::NEUMANN_CONSTANT_SOURCES_AND_VORTEX);
//	nde::Vector<double> F3 = body2D.getForceCoeffs();
//	cout << "Force Neumann = (" << F3(0) << "," << F3(1) << ")" << endl;


//	cout << "Angle of attack = 0.0873" << endl;

//	body2D.calcPotentialFlow(5.0 * M_PI / 180.0,
//									 nde::DIRICHLET_CONSTANT_DOUBLETS);
//	nde::Vector<double> F1_1 = body2D.getForceCoeffs();
//	cout << "Force Dirichlet 1 = (" << F1_1(0) << ", " << F1_1(1) << ")" << endl;

//	body2D.calcPotentialFlow(5.0 * M_PI / 180.0,
//									 nde::DIRICHLET_CONSTANT_SOURCES_AND_DOUBLETS);
//	nde::Vector<double> F2_1 = body2D.getForceCoeffs();
//	cout << "Force Dirichlet 2 = (" << F2_1(0) << "," << F2_1(1) << ")" << endl;

//	body2D.calcPotentialFlow(5.0 * M_PI / 180.0,
//									 nde::NEUMANN_CONSTANT_SOURCES_AND_VORTEX);
//	nde::Vector<double> F3_1 = body2D.getForceCoeffs();
//	cout << "Force Neumann = (" << F3_1(0) << "," << F3_1(1) << ")" << endl;

//	}


//	{

//	cout << "/******************************************************/" << endl;
//	cout << "/* File Input airfoil test" << endl;
//	cout << "/******************************************************/" << endl;

//	nde::Airfoil airfoil("ex_airfoil1.dat", 0.1);
//	airfoil.writeCoordsToFile("ex_airfoil1_out.dat",50);

//	nde::Vector<nde::Panel2D> panels = airfoil.getPanels(0.01);

//	cout << "Angle of attack = 0.0" << endl;
//	nde::AerodynamicBody2D body2D(1.0, panels);

//	body2D.calcPotentialFlow(0.0, nde::DIRICHLET_CONSTANT_DOUBLETS);
//	nde::Vector<double> F1 = body2D.getForceCoeffs();
//	cout << "Force Dirichlet 1 = (" << F1(0) << ", " << F1(1) << ")" << endl;

//	body2D.calcPotentialFlow(0.0, nde::DIRICHLET_CONSTANT_SOURCES_AND_DOUBLETS);
//	nde::Vector<double> F2 = body2D.getForceCoeffs();
//	cout << "Force Dirichlet 2 = (" << F2(0) << "," << F2(1) << ")" << endl;

//	body2D.calcPotentialFlow(0.0, nde::NEUMANN_CONSTANT_SOURCES_AND_VORTEX);
//	nde::Vector<double> F3 = body2D.getForceCoeffs();
//	cout << "Force Neumann = (" << F3(0) << "," << F3(1) << ")" << endl;

//	}


//	{

//	cout << "/******************************************************/" << endl;
//	cout << "/* Interpolation test" << endl;
//	cout << "/******************************************************/" << endl;

//	nde::Vector<double> x(5);
//	nde::Vector<double> y(5);

//	for (int i = 0; i < 5; ++i) {
//		x(i) = double(i);
//		y(i) = double(2 * i + 1);
//	}

//	nde::Interpolator1D lin_interp(x, y, nde::LINEAR);

//	cout << "interpolation at 0.5 = " << lin_interp(0.5)
//		  << ", slope = " << lin_interp(1, 0.5) << endl;
//	cout << "interpolation at 2.5 = " << lin_interp(2.5)
//		  << ", slope = " << lin_interp(1, 2.5) << endl;
//	cout << "interpolation at 3.5 = " << lin_interp(3.5)
//		  << ", slope = " << lin_interp(1, 3.5) << endl;

//	nde::Interpolator1D spln_interp(x, y, nde::SPLINE_MONOTONE);

//	cout << "interpolation at 0.5 = " << spln_interp(0.5)
//		  << ", slope = " << spln_interp(1, 0.5) << endl;
//	cout << "interpolation at 2.5 = " << spln_interp(2.5)
//		  << ", slope = " << spln_interp(1, 2.5) << endl;
//	cout << "interpolation at 3.5 = " << spln_interp(3.5)
//		  << ", slope = " << spln_interp(1, 3.5) << endl;

//	}



	return 0;

}
