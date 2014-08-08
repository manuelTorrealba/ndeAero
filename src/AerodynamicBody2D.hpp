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
                double chord_in,
                const Vector<Panel2D>& panels_in,
                double angle_attack_in);
        void calcPotentialFlow();
        double getPotential(const Vector<double>& x) const;
        Vector<double> getSpeed(const Vector<double>& x) const;
        Vector<double> AerodynamicBody2D::calcForceCoefficients() const;

    private:
        double chord;
        double angle_attack;
        Vector<Panel2D> panels;
        Vector<double> incident_flow;
        Vector<double> sources;
        Vector<double> doublets;
        Vector<double> vortexes;
        double wake;
        
        double incidentFlowPotential(Vector<double> x) const;
        
    };

}


#endif	/* AERODYNAMICBODY2D_HPP */

