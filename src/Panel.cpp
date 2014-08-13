/* 
 * File:   Panel.hpp
 * Author: kenobi
 *
 * Created on July 21, 2014, 9:34 AM
 */

#include"Panel.hpp"

namespace nde {

    Panel2D::Panel2D(
            const Vector<double>& start_point_in,
            const Vector<double>& end_point_in)
    : start_point(start_point_in), end_point(end_point_in) {

        setPanel2D();

    }

    void Panel2D::setPoints(const Vector<double>& start_point_in,
            const Vector<double>& end_point_in) {

        start_point = start_point_in;
        end_point = end_point_in;
        setPanel2D();

    }

    void Panel2D::setPanel2D() {

        dl = end_point - start_point;
        length = dl.norm();
        tangent = dl / length;
        mid_point = (start_point + end_point) * 0.5;
        normal.resize(2);
        normal(0) = -tangent(1);
        normal(1) = tangent(0);

        double theta = atan2(tangent(1), tangent(0));
        B.resize(2, 2);
        B(0, 0) = cos(theta);
        B(0, 1) = sin(theta);
        B(1, 0) = -sin(theta);
        B(1, 1) = cos(theta);

        start_point_local = B*start_point;
        end_point_local = B*end_point;

        control_point_in = mid_point - normal * (length * 1.0e-4);
        control_point_out = mid_point + normal * (length * 1.0e-4);

    }

    Vector<double> Panel2D::getStartPoint() const {
        return start_point;
    }

    Vector<double> Panel2D::getEndPoint() const {
        return end_point;
    }

    double Panel2D::getDx() const {
        return dl(0);
    }

    double Panel2D::getDz() const {
        return dl(1);
    }

    Vector<double> Panel2D::getDl() const {
        return dl;
    }

    double Panel2D::getLength() const {
        return length;
    }

    Vector<double> Panel2D::getMidPoint() const {
        return mid_point;
    }

    Vector<double> Panel2D::getControlPointIn() const {
        return control_point_in;
    }

    Vector<double> Panel2D::getControlPointOut() const {
        return control_point_out;
    }

    Vector<double> Panel2D::getNormal() const {
        return normal;
    }

    Vector<double> Panel2D::getTangent() const {
        return tangent;
    }

    double Panel2D::calcConstantDoubletPotencial(const Vector<double>& x) const {

        Vector<double> xl = B*x;

        double theta1 = atan2(xl(1) - start_point_local(1), xl(0) - start_point_local(0));
        double theta2 = atan2(xl(1) - end_point_local(1), xl(0) - end_point_local(0));

        return -(theta2 - theta1) / (2.0 * M_PI);

    }

    double Panel2D::calcConstantSourcePotencial(const Vector<double>& x) const {

        double r1 = (x - start_point).norm();
        double r2 = (x - end_point).norm();
        Vector<double> xl = B*x;
        double theta1 = atan2(xl(1) - start_point_local(1), xl(0) - start_point_local(0));
        double theta2 = atan2(xl(1) - end_point_local(1), xl(0) - end_point_local(0));

        double tol = 1.0e-6 * length;

        double t1;
        if (r1 > tol)
            t1 = 2.0 * r1 * cos(theta1) * log(r1);
        else
            t1 = 0.0;

        double t2;
        if (r2 > tol)
            t2 = -2.0 * r2 * cos(theta2) * log(r2);
        else
            t2 = 0.0;

        double t3 = 2.0 * r1 * sin(theta1)*(theta2 - theta1);

        return (t1 + t2 + t3) / (4.0 * M_PI);

    }

    Vector<double> Panel2D::calcConstantDoubletSpeed(const Vector<double>& x) const {

        Vector<double> u(2);
        Vector<double> r1 = x - start_point;
        Vector<double> r2 = x - end_point;
        double r1sqr = r1.norm() * r1.norm();
        double r2sqr = r2.norm() * r2.norm();

        u(0) = 1.0 / (2.0 * M_PI)* (r2(1) / r2sqr - r1(1) / r1sqr);
        u(1) = -1.0 / (2.0 * M_PI)* (r2(0) / r2sqr - r1(0) / r1sqr);

        return u;

    }

    Vector<double> Panel2D::calcConstantSourceSpeed(const Vector<double>& x) const {

        Vector<double> xl = B*x;
        double theta1 = atan2(xl(1) - start_point_local(1), xl(0) - start_point_local(0));
        double theta2 = atan2(xl(1) - end_point_local(1), xl(0) - end_point_local(0));
        Vector<double> r1 = x - start_point;
        Vector<double> r2 = x - end_point;

        Vector<double> u(2);
        u(0) = 1.0 / (2.0 * M_PI) * log(r1.norm() / r2.norm());
        u(1) = 1.0 / (2.0 * M_PI) * (theta2 - theta1);

        return B.solve(u); //return in global coordinates

    }

    Vector<double> Panel2D::calcConstantVortexSpeed(const Vector<double>& x) const {

        Vector<double> xl = B*x;
        double theta1 = atan2(xl(1) - start_point_local(1), xl(0) - start_point_local(0));
        double theta2 = atan2(xl(1) - end_point_local(1), xl(0) - end_point_local(0));
        Vector<double> r1 = x - start_point;
        Vector<double> r2 = x - end_point;

        Vector<double> u(2);
        u(0) = 1.0 / (2.0 * M_PI) * (theta2 - theta1);
        u(1) = 1.0 / (2.0 * M_PI) * log(r2.norm() / r1.norm());

        return B.solve(u); //return in global coordinates

    }


}
