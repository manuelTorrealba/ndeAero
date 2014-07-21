/* 
 * File:   PotentialFlowElements.hpp
 * Author: kenobi
 *
 * Created on July 17, 2014, 11:09 AM
 */

#include"PotentialFlowElements.hpp"

namespace nde {

    namespace potential_flow {

        PointLocationPanel2D::PointLocationPanel2D(const Vector<double>& x1,
                const Vector<double>& x2,
                const Vector<double>& x) {

            Vector<double> d1 = x - x1;
            Vector<double> d2 = x - x2;
            r1 = d1.norm();
            r2 = d2.norm();
            length = (x1 - x2).norm();

            double a1 = acos(-(r2 * r2 - r1 * r1 - length * length) / (2 * r1 * length));
            double a2 = M_PI - acos(-(r1 * r1 - r2 * r2 - length * length) / (2 * r2 * length));

            double p = (x2(1) - x1(1)) / (x2(0) - x1(0));
            double m = (x1(1) - x1(0) * p) / (x(1) - x(0) * p);

            if (m < 1.0) {
                angle1 = a1;
                angle2 = a2;
            } else {
                angle1 = 2.0 * M_PI - a1;
                angle2 = 2.0 * M_PI - a2;
            }

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

        double ConstantDoublet2D_potential(Vector<double> x1, Vector<double> x2, Vector<double> x) {
            PointLocationPanel2D Loc(x1, x2, x);
            return -(Loc.getAngle2() - Loc.getAngle1()) / (2.0 * M_PI);
        }

        double PointVortex2D_potential(Vector<double> x0, Vector<double> x) {
            return -atan((x(1) - x0(1)) / (x(0) - x0(0))) / (2.0 * M_PI);
        }

    }

}

