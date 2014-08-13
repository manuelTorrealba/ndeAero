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

    enum PanelMethodType {
        DIRICHLET_CONSTANT_DOUBLETS = 1,
        NEUMANN_CONSTANT_SOURCES_AND_VORTEX = 2
    };
    
    class AerodynamicBody2D {
    public:
        AerodynamicBody2D(
                double chord_in,
                const Vector<Panel2D>& panels_in,
                double angle_attack_in);
        void calcPotentialFlow(PanelMethodType panel_method_type);
        Vector<double> getForceCoeffs() const;
        
    private:
        double chord;
        double angle_attack;
        Vector<Panel2D> panels;
        Vector<double> incident_flow;

        Vector<double> x;
        Vector<double> v;
        Vector<double> cp;
        Vector<double> F;

    };

}


#endif	/* AERODYNAMICBODY2D_HPP */

