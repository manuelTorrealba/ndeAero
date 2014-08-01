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

        control_point_in = mid_point - normal * (length * 1.0e-4);
        control_point_out = mid_point + normal * (length * 1.0e-4);
    }

    const Vector<double>& Panel2D::getStartPoint() const {
        return start_point;
    }

    const Vector<double>& Panel2D::getEndPoint() const {
        return end_point;
    }

    const double Panel2D::getDx() const {
        return dl(0);
    }

    const double Panel2D::getDz() const {
        return dl(1);
    }

    const Vector<double>& Panel2D::getDl() const {
        return dl;
    }

    const double Panel2D::getLength() const {
        return length;
    }

    const Vector<double>& Panel2D::getMidPoint() const {
        return mid_point;
    }

    const Vector<double>& Panel2D::getControlPointIn() const {
        return control_point_in;
    }

    const Vector<double>& Panel2D::getControlPointOut() const {
        return control_point_out;
    }

    const Vector<double>& Panel2D::getNormal() const {
        return normal;
    }

    const Vector<double>& Panel2D::getTangent() const {
        return tangent;
    }

}
