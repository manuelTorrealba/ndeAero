#include"Vector.hpp"
#include"Matrix.hpp"
#include"PotentialFlowElements.hpp"
#include"AerodynamicBody2D.hpp"
#include"Airfoil.hpp"

#include<iostream>


/* nde: numerical design engineering */


using namespace std;

int main(int narg, char** arg) {

    {

        nde::Vector<double> v(5);

        for (int i = 0; i < 5; ++i) {
            v(i) = double(i);
            cout << v(i) << endl;
        }

        cout << v * v << endl;

        v.resize(2);
        v(0) = 10;
        v(1) = 20;
        cout << v(0) << "," << v(1) << ",size=" << v.size() << endl;

        nde::Matrix<double> A(2, 2);
        A(0, 0) = 1.0;
        A(0, 1) = 2.0;
        A(1, 0) = 3.0;
        A(1, 1) = 4.0;

        nde::Matrix<double> B(2, 2);
        B(0, 0) = 2.0;
        B(0, 1) = 1.0;
        B(1, 0) = 1.0;
        B(1, 1) = 1.0;

        cout << "1" << endl;
        /* vector by matrix multiplication */
        cout << "number of rows=" << A.numrows() << endl;
        cout << "number of cols=" << A.numcols() << endl;

        cout << A(0, 0) << "," << A(0, 1) << endl;
        cout << A(1, 0) << "," << A(1, 1) << endl;

        nde::Vector<double> u = A.row(0);

        cout << u(0) << "," << u(1) << endl;

        nde::Vector<double> w = A*v;
        cout << "2" << endl;

        cout << "vector by matrix:" << endl;
        cout << "Result=(" << w(0) << "," << w(1) << ")" << endl;
        cout << "Expected=(50,110)" << endl;

        nde::Matrix<double> Mat(3, 3);
        Mat(0, 0) = 3.0;
        Mat(0, 1) = 1.0;
        Mat(0, 2) = -1.0;
        Mat(1, 0) = 1.0;
        Mat(1, 1) = 0.5;
        Mat(1, 2) = 1.0;
        Mat(2, 0) = 2.0;
        Mat(2, 1) = -5.0;
        Mat(2, 2) = 2.0;

        nde::Vector<double> b(3);
        b(0) = 2.0;
        b(1) = 5.0;
        b(2) = -2.0;

        nde::Vector<double> sol = Mat.solve(b);

        cout << "Solution of the linear system =( "
                << sol(0) << "," << sol(1) << "," << sol(2) << ")" << endl;


    }


    /* Constant source 2D potential test */
    {
        nde::Vector<double> x(2);
        x(0) = -0.5;
        x(1) = 1.0;

        nde::Vector<double> x1(2);
        x1(0) = -0.5;
        x1(1) = 0.0;

        nde::Vector<double> x2(2);
        x2(0) = 0.5;
        x2(1) = 0.0;

        double phi_source2D = nde::potential_flow::ConstantSource2D_potential(x1, x2, x);

        cout << "Constant Source 2D potential" << endl;
        cout << "Result = " << phi_source2D << endl;
        cout << "Expected = " << (log(2.0) + M_PI / 2.0) / (4.0 * M_PI) << endl;

    }


    /* Constant doublet 2D potential test */
    {

        nde::Vector<double> x(2);
        x(0) = 1.0;
        x(1) = 0.0;

        nde::Vector<double> x1(2);
        x1(0) = 0.0;
        x1(1) = 0.0;

        nde::Vector<double> x2(2);
        x2(0) = 0.0;
        x2(1) = 1.0;

        double phi_doublet2D = nde::potential_flow::ConstantDoublet2D_potential(x1, x2, x);

        cout << "Constant Doublet 2D potential" << endl;
        cout << "Result = " << phi_doublet2D << endl;
        cout << "Expected = " << 0.125 << endl;

    }

    /* Vortex potential test */
    {

        nde::Vector<double> x(2);
        x(0) = 2.5;
        x(1) = 2.5;

        nde::Vector<double> x1(2);
        x1(0) = 1.0;
        x1(1) = 1.0;

        double phi_vortex2D = nde::potential_flow::PointVortex2D_potential(x1, x);

        cout << "Constant Vortex 2D potential" << endl;
        cout << "Result = " << phi_vortex2D << endl;
        cout << "Expected = " << -0.125 << endl;

    }


    /* Fluid Elements test */
    {

        nde::Vector<double> x0(2);
        x0.fill(0.0);

        nde::Vector<double> x1(2);
        x1(0) = 0.05 * cos(210. * M_PI / 180.);
        x1(1) = 0.05 * sin(210. * M_PI / 180.);

        nde::Vector<double> x2(2);

        x2(0) = 0.05 * cos(30. * M_PI / 180.);
        x2(1) = 0.05 * sin(30. * M_PI / 180.);

        for (int i = 0; i < 8; ++i) {
            double angle = 30.0 + double(i) * 360. / 8.;
            nde::Vector<double> x(2);
            x(0) = 1.0 * cos(angle * M_PI / 180.);
            x(1) = 1.0 * sin(angle * M_PI / 180.);
            cout << "Potential Source, angle " << angle << " = " <<
                    nde::potential_flow::ConstantSource2D_potential(x1, x2, x) << endl;
        }

        for (int i = 0; i < 8; ++i) {
            double angle = 30.0 + double(i) * 360. / 8.;
            nde::Vector<double> x(2);
            x(0) = 1.0 * cos(angle * M_PI / 180.);
            x(1) = 1.0 * sin(angle * M_PI / 180.);
            nde::Vector<double> v = nde::potential_flow::ConstantSource2D_speed(x1, x2, x);
            cout << "Speed Source, angle " << angle << " = " << v(0) << "," << v(1) << endl;
        }

        for (int i = 0; i < 8; ++i) {
            double angle = 30.0 + double(i) * 360. / 8.;
            nde::Vector<double> x(2);
            x(0) = 1.0 * cos(angle * M_PI / 180.);
            x(1) = 1.0 * sin(angle * M_PI / 180.);
            cout << "Potential Doublet, angle " << angle << " = " <<
                    nde::potential_flow::ConstantDoublet2D_potential(x1, x2, x) << endl;
        }

        for (int i = 0; i < 8; ++i) {
            double angle = 30.0 + double(i) * 360. / 8.;
            nde::Vector<double> x(2);
            x(0) = 1.0 * cos(angle * M_PI / 180.);
            x(1) = 1.0 * sin(angle * M_PI / 180.);
            nde::Vector<double> v = nde::potential_flow::ConstantDoublet2D_speed(x1, x2, x);
            cout << "Speed Doublet, angle " << angle << " = " << v(0) << "," << v(1) << endl;
        }

        for (int i = 0; i < 8; ++i) {
            double angle = 30.0 + double(i) * 360. / 8.;
            nde::Vector<double> x(2);
            x(0) = 1.0 * cos(angle * M_PI / 180.);
            x(1) = 1.0 * sin(angle * M_PI / 180.);
            cout << "Potential Doublet, angle " << angle << " = " <<
                    nde::potential_flow::PointVortex2D_potential(x2, x) -
                    nde::potential_flow::PointVortex2D_potential(x1, x) << endl;
        }

        for (int i = 0; i < 8; ++i) {
            double angle = 30.0 + double(i) * 360. / 8.;
            nde::Vector<double> x(2);
            x(0) = 1.0 * cos(angle * M_PI / 180.);
            x(1) = 1.0 * sin(angle * M_PI / 180.);
            nde::Vector<double> v = nde::potential_flow::PointVortex2D_speed(x2, x) -
                    nde::potential_flow::PointVortex2D_speed(x1, x);
            cout << "Speed Doublet, angle " << angle << " = " << v(0) << "," << v(1) << endl;
        }

        for (int i = 0; i < 8; ++i) {
            double angle = 30.0 + double(i) * 360. / 8.;
            nde::Vector<double> x(2);
            x(0) = 1.0 * cos(angle * M_PI / 180.);
            x(1) = 1.0 * sin(angle * M_PI / 180.);
            cout << "Potential Vortex, angle " << angle << " = " <<
                    nde::potential_flow::PointVortex2D_potential(x0, x) << endl;
        }

        for (int i = 0; i < 8; ++i) {
            double angle = 30.0 + double(i) * 360. / 8.;
            nde::Vector<double> x(2);
            x(0) = 1.0 * cos(angle * M_PI / 180.);
            x(1) = 1.0 * sin(angle * M_PI / 180.);
            nde::Vector<double> v = nde::potential_flow::PointVortex2D_speed(x0, x);
            cout << "Speed Vortex, angle " << angle << " = " << v(0) << "," << v(1) << endl;
        }

    }


    /* Aerodynamic body test */
    {
        int num_panels = 10;
        nde::Vector<nde::Panel2D> panels(2 * num_panels);
        for (int i = 0; i < num_panels; ++i) {
            nde::Vector<double> x1(2);
            nde::Vector<double> x2(2);
            x1(0) = 1.0 * i / double(num_panels);
            x2(0) = 1.0 * (i + 1) / double(num_panels);
            if (x1(0) <= 0.5)
                x1(1) = 0.1 * x1(0) / 0.5;
            else
                x1(1) = 0.2 - 0.1 * x1(0) / 0.5;

            if (x2(0) <= 0.5)
                x2(1) = 0.1 * x2(0) / 0.5;
            else
                x2(1) = 0.2 - 0.1 * x2(0) / 0.5;

            panels(i).setPoints(x1, x2);

        }

        for (int i = 0; i < num_panels; ++i) {
            nde::Vector<double> x1(2);
            nde::Vector<double> x2(2);
            x1(0) = 1.0 - 1.0 * i / double(num_panels);
            x2(0) = 1.0 - 1.0 * (i + 1) / double(num_panels);
            if (x1(0) <= 0.5)
                x1(1) = 0.0; //-0.1 * x1(0) / 0.5;
            else
                x1(1) = 0.0; //-0.2 + 0.1 * x1(0) / 0.5;

            if (x2(0) <= 0.5)
                x2(1) = 0.0; //-0.1 * x2(0) / 0.5;
            else
                x2(1) = 0.0; //-0.2 + 0.1 * x2(0) / 0.5;

            panels(i + num_panels).setPoints(x1, x2);

        }

        nde::Vector<double> wake_coordinates(2);
        wake_coordinates(0) = 1.0;
        wake_coordinates(1) = 0.0;
        double air_speed = 1.0;
        double angle_attack = 0.0;

        nde::AerodynamicBody2D body2D(
                panels, wake_coordinates, air_speed, angle_attack);

        body2D.calcPotentialFlow();

        nde::Vector<nde::Vector<double> > p_out(4);

        p_out(0).resize(2);
        p_out(0)(0) = -1.0;
        p_out(0)(1) = 0.0;

        p_out(1).resize(2);
        p_out(1)(0) = 0.5;
        p_out(1)(1) = 1.0;

        p_out(2).resize(2);
        p_out(2)(0) = 2.0;
        p_out(2)(1) = 1.0;

        p_out(3).resize(2);
        p_out(3)(0) = 0.5;
        p_out(3)(1) = -1.0;

        for (int i = 0; i < panels.size(); ++i) {
            nde::Vector<double> control_point =
                    panels(i).getMidPoint() - panels(i).getNormal()*(panels(i).getLength()* 1.0e-3);
            double phi_control = body2D.getPotential(control_point);
            cout << "Control point (" << i + 1 << ") = " << control_point(0)
                    << "," << control_point(1) << endl;
            cout << "Potential control point (" << i + 1 << ") = " << phi_control << endl;
            nde::Vector<double> u_control = body2D.getSpeed(
                    panels(i).getMidPoint() + panels(i).getNormal()*(panels(i).getLength()* 1.0e-3));
            cout << "Speed control point (" << i + 1 << ") = " << u_control(0)
                    << "," << u_control(1) << endl;
        }

        for (int i = 0; i < p_out.size(); ++i) {
            double phi_out = body2D.getPotential(p_out(i));
            cout << "Potential point out (" << i + 1 << ") = " << phi_out << endl;
        }

    }



    {

        nde::NacaAirfoil naca_airfoil(1.0, 0, 0, 2, 0);

        nde::Airfoil airfoil(naca_airfoil);

        nde::Vector<nde::Panel2D> panels = airfoil.getPanels(0.2);

        for (int i = 0; i < panels.size(); ++i) {
            nde::Vector<double> x = panels(i).getMidPoint();
            cout << "Panel mid point " << i << " = " << x(0) << "," << x(1) << endl;
        }

    }

    return 0;

}
