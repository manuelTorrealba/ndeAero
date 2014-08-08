/* 
 * File:   PotentialFlowElements.hpp
 * Author: kenobi
 *
 * Created on July 17, 2014, 11:09 AM
 */

#include <math.h>

#include"PotentialFlowElements.hpp"

namespace nde {

    namespace potential_flow {

        PointLocationPanel2D::PointLocationPanel2D(const Vector<double>& x1,
                const Vector<double>& x2,
                const Vector<double>& x) {

            Vector<double> dx = x2 - x1;
            theta = atan2(dx(1), dx(0));
            length = dx.norm();

            Vector<double> d1 = x - x1;
            Vector<double> d1p(2);
            d1p(0) = d1(0) * cos(theta) + d1(1) * sin(theta);
            d1p(1) = -d1(0) * sin(theta) + d1(1) * cos(theta);
            r1 = d1.norm();
            angle1 = atan2(d1p(1), d1p(0));

            Vector<double> d2 = x - x2;
            Vector<double> d2p(2);
            d2p(0) = d2(0) * cos(theta) + d2(1) * sin(theta);
            d2p(1) = -d2(0) * sin(theta) + d2(1) * cos(theta);
            r2 = d2.norm();
            angle2 = atan2(d2p(1), d2p(0));

        }

        double PointLocationPanel2D::getR1() const {
            return r1;
        }

        double PointLocationPanel2D::getR2() const {
            return r2;
        }

        double PointLocationPanel2D::getLength() const {
            return length;
        }

        double PointLocationPanel2D::getAngle1() const {
            return angle1;
        }

        double PointLocationPanel2D::getAngle2() const {
            return angle2;
        }

        Vector<double> PointLocationPanel2D::fromLocalToGlobalCoordinates(Vector<double> u_l) const {

            Vector<double> u_g(2);

            u_g(0) = u_l(0) * cos(theta) - u_l(1) * sin(theta);
            u_g(1) = u_l(0) * sin(theta) + u_l(1) * cos(theta);

            return u_g;

        }

        double ConstantSource2D_potential(Vector<double> x1, Vector<double> x2, Vector<double> x) {

            PointLocationPanel2D Loc(x1, x2, x);
            double length = Loc.getLength();
            double r1 = Loc.getR1();
            double r2 = Loc.getR2();
            double angle1 = Loc.getAngle1();
            double angle2 = Loc.getAngle2();

            double tol = 1.0e-6 * length;

            double t1;
            if (r1 > tol)
                t1 = 2.0 * r1 * cos(angle1) * log(r1);
            else
                t1 = 0.0;

            double t2;
            if (r2 > tol)
                t2 = -2.0 * r2 * cos(angle2) * log(r2);
            else
                t2 = 0.0;

            double t3 = 2.0 * r1 * sin(angle1)*(angle2 - angle1);

            return (t1 + t2 + t3) / (4.0 * M_PI);

        }

        Vector<double> ConstantSource2D_speed(Vector<double> x1, Vector<double> x2, Vector<double> x) {

            PointLocationPanel2D Loc(x1, x2, x);

            double r1 = Loc.getR1();
            double r2 = Loc.getR2();
            double angle1 = Loc.getAngle1();
            double angle2 = Loc.getAngle2();

            Vector<double> u(2);
            u(0) = 1.0 / (2.0 * M_PI) * log(r1 / r2);
            u(1) = 1.0 / (2.0 * M_PI) * (angle2 - angle1);

            return Loc.fromLocalToGlobalCoordinates(u);

        }

        double ConstantDoublet2D_potential(Vector<double> x1, Vector<double> x2, Vector<double> x) {
            PointLocationPanel2D Loc(x1, x2, x);
            return -(Loc.getAngle2() - Loc.getAngle1()) / (2.0 * M_PI);
        }

        Vector<double> ConstantDoublet2D_speed(Vector<double> x1, Vector<double> x2, Vector<double> x) {

            PointLocationPanel2D Loc(x1, x2, x);

            double r1 = Loc.getR1();
            double r2 = Loc.getR2();
            double angle1 = Loc.getAngle1();
            double angle2 = Loc.getAngle2();

            Vector<double> u(2);
            u(0) = -1.0 / (2.0 * M_PI)*(sin(angle1) / r1 - sin(angle2) / r2);
            u(1) = 1.0 / (2.0 * M_PI)*(cos(angle1) / r1 - cos(angle2) / r2);

            return Loc.fromLocalToGlobalCoordinates(u);

        }

        double ConstantVortex2D_potential(Vector<double> x1, Vector<double> x2, Vector<double> x) {
            return 0.0;
        }

        Vector<double> ConstantVortex2D_speed(Vector<double> x1, Vector<double> x2, Vector<double> x) {

            PointLocationPanel2D Loc(x1, x2, x);

            double r1 = Loc.getR1();
            double r2 = Loc.getR2();
            double angle1 = Loc.getAngle1();
            double angle2 = Loc.getAngle2();

            Vector<double> u(2);
            u(0) = 1.0 / (2.0 * M_PI) * (angle2 - angle1);
            u(1) = 1.0 / (2.0 * M_PI) * log(r2 / r1);

            return Loc.fromLocalToGlobalCoordinates(u);

        }

        double PointVortex2D_potential(Vector<double> x0, Vector<double> x) {
            return -atan((x(1) - x0(1)) / (x(0) - x0(0))) / (2.0 * M_PI);
        }

        Vector<double> PointVortex2D_speed(Vector<double> x0, Vector<double> x) {

            Vector<double> u(2);
            double dz = x(1) - x0(1);
            double dx = x(0) - x0(0);

            u(0) = 1.0 / (2 * M_PI) * dz / (dx * dx + dz * dz);
            u(1) = -1.0 / (2 * M_PI) * dx / (dx * dx + dz * dz);

            return u;

        }

    }

}

