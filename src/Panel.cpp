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
        dl = end_point - start_point;
        length = dl.norm();
        mid_point = (start_point + end_point) * 0.5;
        Vector<double> n(2);
        n(0) = -dl(1);
        n(1) = dl(0);
        normal = n / n.norm();
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

    const Vector<double>& Panel2D::getNormal() const {
        return normal;
    }

}
