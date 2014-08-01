/* 
 * File:   Panel.hpp
 * Author: kenobi
 *
 * Created on July 21, 2014, 9:34 AM
 */

#ifndef PANEL_HPP
#define	PANEL_HPP

#include"Vector.hpp"

namespace nde {

    class Panel2D {
    public:
        Panel2D(const Vector<double>& start_point_in,
                const Vector<double>& end_point_in);

        void setPoints(const Vector<double>& start_point_in,
                const Vector<double>& end_point_in);

        const Vector<double>& getStartPoint() const;
        const Vector<double>& getEndPoint() const;
        const double getDx() const;
        const double getDz() const;
        const Vector<double>& getDl() const;
        const double getLength() const;
        const Vector<double>& getMidPoint() const;
        const Vector<double>& getControlPointIn() const;
        const Vector<double>& getControlPointOut() const;
        const Vector<double>& getNormal() const;
        const Vector<double>& getTangent() const;

    private:
        Vector<double> start_point;
        Vector<double> end_point;
        Vector<double> dl;
        double length;
        Vector<double> mid_point;
        Vector<double> control_point_in;
        Vector<double> control_point_out;
        Vector<double> normal;
        Vector<double> tangent;
        void setPanel2D();
    };

}


#endif	/* PANEL_HPP */

