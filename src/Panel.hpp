/* 
 * File:   Panel.hpp
 * Author: kenobi
 *
 * Created on July 21, 2014, 9:34 AM
 */

#ifndef PANEL_HPP
#define	PANEL_HPP

#include"Vector.hpp"
#include"Matrix.hpp"

namespace nde {

    class Panel2D {
    public:
        Panel2D(const Vector<double>& start_point_in,
                const Vector<double>& end_point_in);

        void setPoints(const Vector<double>& start_point_in,
                const Vector<double>& end_point_in);

        Vector<double> getStartPoint() const;
        Vector<double> getEndPoint() const;
        double getDx() const;
        double getDz() const;
        Vector<double> getDl() const;
        double getLength() const;
        Vector<double> getMidPoint() const;
        Vector<double> getControlPointIn() const;
        Vector<double> getControlPointOut() const;
        Vector<double> getNormal() const;
        Vector<double> getTangent() const;

        double calcConstantDoubletPotencial(const Vector<double>& x) const;
        
        Vector<double> calcConstantSourceSpeed(const Vector<double>& x) const;
        Vector<double> calcConstantDoubletSpeed(const Vector<double>& x) const;
        Vector<double> Panel2D::calcConstantVortexSpeed(const Vector<double>& x) const;
    
    private:
        Vector<double> start_point;
        Vector<double> end_point;
        Vector<double> start_point_local;
        Vector<double> end_point_local;
        Vector<double> dl;
        double length;
        Vector<double> mid_point;
        Vector<double> control_point_in;
        Vector<double> control_point_out;
        Vector<double> normal;
        Vector<double> tangent;
        Matrix<double> B; //matrix to change from global to local coordinates
        void setPanel2D();

    };

}


#endif	/* PANEL_HPP */

