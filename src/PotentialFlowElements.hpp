/* 
 * File:   PotentialFlowElements.hpp
 * Author: kenobi
 *
 * Created on July 17, 2014, 11:09 AM
 */

#ifndef POTENTIALFLOWELEMENTS_HPP
#define	POTENTIALFLOWELEMENTS_HPP

#include"Vector.hpp"

namespace nde {

    namespace potential_flow {

        class PointLocationPanel2D {
        public:

            PointLocationPanel2D(const Vector<double>& x1,
                                 const Vector<double>& x2,
                                 const Vector<double>& x);

            double getR1() const;
            double getR2() const;
            double getLength() const;
            double getAngle1() const;
            double getAngle2() const;
            Vector<double> fromLocalToGlobalCoordinates(Vector<double> u_l) const;
            
        private:
            double r1;
            double r2;
            double length;
            double angle1;
            double angle2;
            double cosAngleLG;
            double sinAngleLG;
        };

        double ConstantSource2D_potential(Vector<double> x1, Vector<double> x2, Vector<double> x);

        Vector<double> ConstantSource2D_speed(Vector<double> x1, Vector<double> x2, Vector<double> x);
        
        double ConstantDoublet2D_potential(Vector<double> x1, Vector<double> x2, Vector<double> x);

        Vector<double> ConstantDoublet2D_speed(Vector<double> x1, Vector<double> x2, Vector<double> x);
        
        double PointVortex2D_potential(Vector<double> x0, Vector<double> x);
        
        Vector<double> PointVortex2D_speed(Vector<double> x0, Vector<double> x);

    }

}

#endif	/* POTENTIALFLOWELEMENTS_HPP */

