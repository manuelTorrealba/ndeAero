#include"Vector.hpp"
#include"Matrix.hpp"
#include"PotentialFlowElements.hpp"
#include"AerodynamicBody2D.hpp"

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

    /* Aerodynamic body test */
    {
        nde::Vector<double> p1(2);
        p1(0) = 0.0;
        p1(1) = 0.0;
        nde::Vector<double> p2(2);
        p2(0) = 0.5;
        p2(1) = 0.1;
        nde::Vector<double> p3(2);
        p3(0) = 1.0;
        p3(1) = 0.0;

        nde::Panel2D panel1(p1, p2);
        nde::Panel2D panel2(p2, p3);
        nde::Panel2D panel3(p3, p1);

        nde::Vector<nde::Panel2D> panels(3);
        panels(0) = panel1;
        panels(1) = panel2;
        panels(2) = panel3;

        nde::Vector<double> wake_coordinates = p3;
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
            double phi_control = body2D.getPotential(
                    panels(i).getMidPoint() - panels(i).getNormal()*1.0e-3);
            cout << "Potential control point (" << i + 1 << ") = " << phi_control << endl;
        }

        for (int i = 0; i < p_out.size(); ++i) {
            double phi_out = body2D.getPotential(p_out(i));
            cout << "Potential point out (" << i + 1 << ") = " << phi_out << endl;
        }
    }

    return 0;

}
