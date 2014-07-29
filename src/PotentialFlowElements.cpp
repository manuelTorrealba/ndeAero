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
            Vector<double> d1 = x - x1;
            Vector<double> d2 = x - x2;
            r1 = d1.norm();
            r2 = d2.norm();
            length = dx.norm();
            Vector<double> dz(2);
            dz(0) = -dx(1);
            dz(1) = dx(0);

            double a1 = acos(-(r2 * r2 - r1 * r1 - length * length) / (2 * r1 * length));
            double a2 = M_PI - acos(-(r1 * r1 - r2 * r2 - length * length) / (2 * r2 * length));

            if (dz * d1 > 0.0) {
                angle1 = a1;
                angle2 = a2;
            } else {
                angle1 = 2.0 * M_PI - a1;
                angle2 = 2.0 * M_PI - a2;
            }

            Vector<double> xG(2); // global x axis.
            xG(0) = 1.0;
            xG(1) = 0.0;
            Vector<double> xL = x2 - x1; // local x axis.
            xL = xL / xL.norm();

            double sgn;
            if (xL(1) >= 0.0)
                sgn = 1.0;
            else
                sgn = -1.0;

            cosAngleLG = xG * xL;
            sinAngleLG = sgn * sqrt(1.0 - cosAngleLG * cosAngleLG);

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

            u_g(0) = u_l(0) * cosAngleLG - u_l(1) * sinAngleLG;
            u_g(1) = u_l(0) * sinAngleLG + u_l(1) * cosAngleLG;

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

        double PointVortex2D_potential(Vector<double> x0, Vector<double> x) {
            double angle = atan2(x(1) - x0(1), x(0) - x0(0));
            if (angle < 0.0) angle = 2 * M_PI + angle;
            return -angle / (2.0 * M_PI);
            //return -atan2((x(1) - x0(1)) / (x(0) - x0(0))) / (2.0 * M_PI);
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

