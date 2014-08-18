/* 
 * File:   Airfoil.hpp
 * Author: root
 *
 * Created on July 23, 2014, 12:51 PM
 */

#ifndef AIRFOIL_HPP
#define	AIRFOIL_HPP

#include"Vector.hpp"
#include"Panel.hpp"

namespace nde {

    struct NacaAirfoil {
    public:
        double chord;
        int naca1;
        int naca2;
        int naca3;
        int naca4;

        NacaAirfoil(double chord_in,
                int naca1_in,
                int naca2_in,
                int naca3_in,
                int naca4_in)
        : chord(chord_in),
        naca1(naca1_in),
        naca2(naca2_in),
        naca3(naca3_in),
        naca4(naca4_in) {
            ;
        }

        double top(double x) const {
            return camber(x) + thickness(x);
        }
        
        double bottom(double x) const {
            return camber(x) - thickness(x);
        }
        
        double camber(double x) const {
            double p = double(naca1) / 100.;
            double m = double(naca2) / 10.;
            double y = 0.0;
            if (x < p * chord)
                y = m * x / (p * p)*(2.0 * p - x / chord);
            else
                y = m * (chord - x) / ((1.0 - p)*(1.0 - p))
                *(1.0 + x / chord - 2.0 * p);
            return y;
        }

        double thickness(double x) const {
            double t = (double(naca3)*10. + double(naca4)) / 100.;
            double h = x / chord;
            double y = t / 0.2 * (0.2969 * sqrt(h) - 0.1260 * h
                    - 0.3516 * h * h + 0.2843 * h * h * h
                    - 0.1036 * h * h * h * h);
            //last coefficient modified from -0.1015 to -0.1036
            //to get zero thickness at the trailing edge while
            //modifying the shape of the airfoil as little
            //as possible.
            return y*chord;
        }

    };

    class Airfoil {
    public:
        Airfoil(const NacaAirfoil& naca_airfoil_in);
        Vector<Panel2D > getPanels(double density) const;
        Vector<double> getTralingEdgeCoordinates() const;
    private:
        NacaAirfoil naca_airfoil;
    };

}

#endif	/* AIRFOIL_HPP */

