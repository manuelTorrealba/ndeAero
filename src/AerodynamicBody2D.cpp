/* 
 * File:   AerodynamicBody2D.hpp
 * Author: kenobi
 *
 * Created on July 21, 2014, 9:29 AM
 */

#include"AerodynamicBody2D.hpp"
#include"PotentialFlowElements.hpp"
#include"Matrix.hpp"

namespace nde {

    AerodynamicBody2D::AerodynamicBody2D(
            const Vector<Panel2D>& panels,
            const Vector<double>& wake_coordinates,
            double air_speed,
            double angle_attack)
    : panels(panels), wake_coordinates(wake_coordinates) {
        incident_flow.resize(2);
        incident_flow(0) = air_speed * cos(angle_attack);
        incident_flow(1) = air_speed * sin(angle_attack);

    }

    void AerodynamicBody2D::calcPotentialFlow() {

        int num_panels = panels.size();

        /* point where the potential = 0 boundary condition
         *  is going to be enforced */
        Vector<Vector<double> > panel_mid_points_in(num_panels);
        for (int i = 0; i < num_panels; ++i)
            panel_mid_points_in(i) = panels(i).getMidPoint() - panels(i).getNormal()*1.0e-3;

        /* calculate the sources intensities */
        sources.resize(num_panels);
        for (int i = 0; i < num_panels; ++i)
            sources(i) = panels(i).getNormal() * incident_flow;

        /* calculate the doublets intensity */

        /* influence coefficients */
        Vector<double> w(num_panels); // wake
        Matrix<double> b(num_panels, num_panels); // sources
        Matrix<double> c(num_panels, num_panels); // doublets

        for (int i = 0; i < num_panels; ++i) {

            // wake
            w(i) = potential_flow::PointVortex2D_potential(
                    wake_coordinates,
                    panel_mid_points_in(i));

            for (int j = 0; j < num_panels; ++j) {

                // sources
                b(i, j) = potential_flow::ConstantSource2D_potential(
                        panels(j).getStartPoint(),
                        panels(j).getEndPoint(),
                        panel_mid_points_in(i));
                // doublets            
                c(i, j) = potential_flow::ConstantDoublet2D_potential(
                        panels(j).getStartPoint(),
                        panels(j).getEndPoint(),
                        panel_mid_points_in(i));
            }

        }


        // system matrix
        Matrix<double> A(num_panels, num_panels);
        for (int i = 0; i < num_panels; ++i) {
            for (int j = 0; j < num_panels; ++j) {
                if (j == 0)
                    A(i, j) = c(i, j) - w(i);
                else if (j == num_panels - 1)
                    A(i, j) = c(i, num_panels - 1) + w(i);
                else
                    A(i, j) = c(i, j);

            }
        }


        Vector<double> rhs = b * (sources * (-1.0));

        // solve for the doublets' intensities
        doublets = A.solve(rhs);

        // wake intensity
        wake = doublets(num_panels - 1) - doublets(0);

    }

    double AerodynamicBody2D::getPotential(const Vector<double>& x) const {

        double phi(0.0);
        for (int i = 0; i < panels.size(); ++i)
            phi += doublets(i) * potential_flow::ConstantDoublet2D_potential(
                panels(i).getStartPoint(),
                panels(i).getEndPoint(),
                x)
            + sources(i) * potential_flow::ConstantSource2D_potential(
                panels(i).getStartPoint(),
                panels(i).getEndPoint(),
                x);

        return phi + wake * potential_flow::PointVortex2D_potential(
                wake_coordinates, x);
        
    }

}
