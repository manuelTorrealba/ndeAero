#include"Vector.hpp"
#include"Matrix.hpp"
#include"Interpolator.hpp"
#include"AerodynamicBody2D.hpp"
#include"Airfoil.hpp"

#include<iostream>

/* nde: numerical design engineering */


using namespace std;

int main(int narg, char** arg) {


    {

        cout << "/******************************************************/" << endl;
        cout << "/* NACA airfoil test" << endl;
        cout << "/******************************************************/" << endl;

        nde::NacaAirfoil naca_airfoil(1.0, 0, 0, 1, 0);
        nde::Airfoil airfoil(naca_airfoil);
        nde::Vector<nde::Panel2D> panels = airfoil.getPanels(0.01);

        nde::AerodynamicBody2D body2D(1.0, panels, 5.0 * M_PI / 180.0);

        body2D.calcPotentialFlow(nde::DIRICHLET_CONSTANT_DOUBLETS);
        nde::Vector<double> F1 = body2D.getForceCoeffs();
        cout << "Force Dirichlet 1 = (" << F1(0) << ", " << F1(1) << ")" << endl;

        body2D.calcPotentialFlow(nde::DIRICHLET_CONSTANT_SOURCES_AND_DOUBLETS);
        nde::Vector<double> F2 = body2D.getForceCoeffs();
        cout << "Force Dirichlet 2 = (" << F2(0) << "," << F2(1) << ")" << endl;

        body2D.calcPotentialFlow(nde::NEUMANN_CONSTANT_SOURCES_AND_VORTEX);
        nde::Vector<double> F3 = body2D.getForceCoeffs();
        cout << "Force Neumann = (" << F3(0) << "," << F3(1) << ")" << endl;

    }

    {

        cout << "/******************************************************/" << endl;
        cout << "/* Interpolation test" << endl;
        cout << "/******************************************************/" << endl;

        nde::Vector<double> x(5);
        nde::Vector<double> y(5);

        for (int i = 0; i < 5; ++i) {
            x(i) = double(i);
            y(i) = double(2 * i + 1);
        }

        nde::Interpolator1D lin_interp(x, y, nde::LINEAR);
        
        cout << "interpolation at 0.5 = " << lin_interp(0.5) << endl;
        cout << "interpolation at 2.5 = " << lin_interp(2.5) << endl;
        cout << "interpolation at 3.5 = " << lin_interp(3.5) << endl;

        nde::Interpolator1D spln_interp(x, y, nde::SPLINE_MONOTONE);
        
        cout << "interpolation at 0.5 = " << spln_interp(0.5) << endl;
        cout << "interpolation at 2.5 = " << spln_interp(2.5) << endl;
        cout << "interpolation at 3.5 = " << spln_interp(3.5) << endl;
        
    }

    return 0;

}
