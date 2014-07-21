/* 
 * File:   AerodynamicBody2D.hpp
 * Author: kenobi
 *
 * Created on July 21, 2014, 9:29 AM
 */

#ifndef AERODYNAMICBODY2D_HPP
#define	AERODYNAMICBODY2D_HPP


#include"Vector.hpp"
#include"Panel.hpp"

namespace nde {

    class AerodynamicBody2D {
    public:
        AerodynamicBody2D(
                const Vector<Panel2D>& panels,
                const Vector<double>& wake_coordinates,
                double air_speed,
                double angle_attack);
        void calcPotentialFlow();
        double getPotential(const Vector<double>& x) const;

    private:
        Vector<Panel2D> panels;
        Vector<double> wake_coordinates;
        Vector<double> incident_flow;
        Vector<double> sources;
        Vector<double> doublets;
        double wake;
        
    };

}


#endif	/* AERODYNAMICBODY2D_HPP */

