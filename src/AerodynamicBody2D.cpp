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
            double chord_in,
            const Vector<Panel2D>& panels_in,
            double angle_attack_in)
    : chord(chord_in), panels(panels_in), angle_attack(angle_attack_in) {
        incident_flow.resize(2);
        incident_flow(0) = cos(angle_attack);
        incident_flow(1) = sin(angle_attack);
    }

    void AerodynamicBody2D::calcPotentialFlow() {

        int num_panels = panels.size();


        /* low speed aerodynamics book program*/
        {
            // find panel angles

            std::cout << "/*****************************************/" << std::endl;
            std::cout << "Book solution" << std::endl;
            std::cout << "/*****************************************/" << std::endl;

            Vector<double> th(num_panels);
            for (int i = 0; i < num_panels; ++i) {
                double dz = panels(i).getEndPoint()(1) - panels(i).getStartPoint()(1);
                double dx = panels(i).getEndPoint()(0) - panels(i).getStartPoint()(0);
                th(i) = atan2(dz, dx);
            }

            // source strengths
            Vector<double> sig(num_panels);

            std::cout << "source strengths (book)" << std::endl;
            for (int i = 0; i < num_panels; ++i) {
                sig(i) = cos(angle_attack) * sin(th(i)) - sin(angle_attack) * cos(th(i));
                //std::cout << sig(i) << std::endl;
            }

            Matrix<double> A(num_panels, num_panels);
            Vector<double> rhs(num_panels);
            Vector<double> w(num_panels);

            /* Influence coefficients */
            for (int i = 0; i < num_panels; ++i) {

                rhs(i) = 0.0;

                for (int j = 0; j < num_panels; ++j) {

                    double xt = panels(i).getMidPoint()(0) - panels(j).getStartPoint()(0);
                    double zt = panels(i).getMidPoint()(1) - panels(j).getStartPoint()(1);

                    double x2t = panels(j).getEndPoint()(0) - panels(j).getStartPoint()(0);
                    double z2t = panels(j).getEndPoint()(1) - panels(j).getStartPoint()(1);

                    double x = xt * cos(th(j)) + zt * sin(th(j));
                    double z = -xt * sin(th(j)) + zt * cos(th(j));
                    double x2 = x2t * cos(th(j)) + z2t * sin(th(j));
                    double z2 = 0.0;

                    double r1 = sqrt(x * x + z * z);
                    double r2 = sqrt((x - x2) * (x - x2) + (z - z2) * (z - z2));
                    double th1 = atan2(z, x);
                    double th2 = atan2(z, x - x2);

                    //doublet influence coefficients
                    if (i == j) {
                        A(i, j) = 0.5;
                        rhs(i) = rhs(i) + sig(j) / 3.14159265 * x * log(r1);
                    } else {
                        A(i, j) = -0.15916 * (th2 - th1);
                        rhs(i) = rhs(i) + sig(j) / 6.28319 *
                                (x * log(r1)-(x - x2) * log(r2) + z * (th2 - th1));
                    }

                }

                w(i) = 0.15916 * atan((panels(i).getMidPoint()(1) - panels(0).getStartPoint()(1)) /
                        (panels(i).getMidPoint()(0) - panels(0).getStartPoint()(0)));

            }


            /* doublets influence coefficients */
            for (int i = 0; i < num_panels; ++i) {
                for (int j = 0; j < num_panels; ++j)
                    std::cout << "A(" << i << "," << j << ")=" << A(i, j) << std::endl;
            }

            /* rhs */
            for (int i = 0; i < num_panels; ++i) {
                std::cout << "rhs(" << i << ")=" << -rhs(i) << std::endl;
            }
            /* wake */
            for (int i = 0; i < num_panels; ++i) {
                std::cout << "wake(" << i << ")=" << w(i) << std::endl;
            }

            Matrix<double> M(num_panels + 1, num_panels + 1);
            Vector<double> B(num_panels + 1);
            M.fill(0.0);
            B.fill(0.0);
            for (int i = 0; i < num_panels; ++i) {
                B(i) = rhs(i);
                M(i, num_panels) = w(i);
                for (int j = 0; j < num_panels; ++j) {
                    M(i, j) = A(i, j);
                }
            }
            M(num_panels, 0) = -1.0;
            M(num_panels, num_panels - 1) = 1.0;
            M(num_panels, num_panels) = -1.0;

            Vector<double> sol = M.solve(B);
            for (int i = 0; i <= num_panels; ++i)
                std::cout << "sol(" << i << ")=" << sol(i) << std::endl;


        }


        std::cout << "/*****************************************/" << std::endl;
        std::cout << "C++ solution" << std::endl;
        std::cout << "/*****************************************/" << std::endl;


        /* calculate the sources intensities */
        sources.resize(num_panels);
        for (int i = 0; i < num_panels; ++i) {
            sources(i) = (-1.0) * (panels(i).getNormal() * incident_flow);
            //std::cout << "sources=" << sources(i) << std::endl;
        }

        /* calculate the doublets intensity */

        /* influence coefficients */
        Vector<double> p(num_panels); // incident flow
        Vector<double> w(num_panels); // wake
        Matrix<double> b(num_panels, num_panels); // sources
        Matrix<double> c(num_panels, num_panels); // doublets


        for (int i = 0; i < num_panels; ++i) {

            // incident flow potential
            p(i) = incidentFlowPotential(panels(i).getControlPointIn());

            // wake
            w(i) = potential_flow::PointVortex2D_potential(
                    panels(0).getStartPoint(),
                    panels(i).getControlPointIn());

            for (int j = 0; j < num_panels; ++j) {

                // sources
                b(i, j) = potential_flow::ConstantSource2D_potential(
                        panels(j).getStartPoint(),
                        panels(j).getEndPoint(),
                        panels(i).getControlPointIn());
                // doublets            
                c(i, j) = potential_flow::ConstantDoublet2D_potential(
                        panels(j).getStartPoint(),
                        panels(j).getEndPoint(),
                        panels(i).getControlPointIn());

            }

        }


        // system matrix
        Matrix<double> A(num_panels, num_panels);
        for (int i = 0; i < num_panels; ++i) {
            for (int j = 0; j < num_panels; ++j) {
                if (j == 0)
                    A(i, j) = c(i, j) - w(i);
                else if (j == num_panels - 1)
                    A(i, j) = c(i, j) + w(i);
                else
                    A(i, j) = c(i, j);

                std::cout << "c(" << i << "," << j << ")=" << c(i, j) << std::endl;

            }
        }


        Vector<double> rhs = b * sources * (-1.0);
        for (int i = 0; i < num_panels; ++i) {
            std::cout << "rhs(" << i << ")=" << rhs(i) << std::endl;
        }

        for (int i = 0; i < num_panels; ++i) {
            std::cout << "wake(" << i << ")=" << w(i) << std::endl;
        }

        // solve for the doublets' intensities
        doublets = A.solve(rhs);

        Vector<double> error = A * doublets - rhs;

        // wake intensity
        wake = doublets(num_panels - 1) - doublets(0);


        //for (int i = 0; i < error.size(); ++i)
        //    std::cout << "Error solution = " << error(i) << std::endl;

        for (int i = 0; i < sources.size(); ++i)
            std::cout << "Sources (" << "i) = " << sources(i) << std::endl;

        for (int i = 0; i < doublets.size(); ++i)
            std::cout << "Doublets (" << "i) = " << doublets(i) << std::endl;

        std::cout << "wake = " << wake << std::endl;

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

        phi = phi + wake * potential_flow::PointVortex2D_potential(
                panels(0).getStartPoint(), x);

        return phi + incidentFlowPotential(x);

    }

    Vector<double> AerodynamicBody2D::getSpeed(const Vector<double>& x) const {

        Vector<double> u(2);
        u.fill(0.0);

        for (int i = 0; i < panels.size(); ++i)
            u = u + potential_flow::ConstantDoublet2D_speed(
                panels(i).getStartPoint(),
                panels(i).getEndPoint(),
                x) * doublets(i)
            + potential_flow::ConstantSource2D_speed(
                panels(i).getStartPoint(),
                panels(i).getEndPoint(),
                x) * sources(i);

        u = u + potential_flow::PointVortex2D_speed(
                panels(0).getStartPoint(), x) * wake;

        return u + incident_flow;

    }

    Vector<double> AerodynamicBody2D::calcForceCoefficients() const {

        Vector<double> cf_global(2);
        cf_global.fill(0.0);

        for (int i = 0; i < panels.size(); ++i) {
            Vector<double> v = getSpeed(panels(i).getControlPointOut());
            cf_global = cf_global + panels(i).getNormal()
                    * ((1.0 - pow(v.norm() , 2))
                    * panels(i).getLength() / chord);
        }

        Vector<double> cf(2);
        cf(0) = cf_global(0) * cos(angle_attack) + cf_global(1) * sin(angle_attack);
        cf(1) = -cf_global(0) * sin(angle_attack) + cf_global(1) * cos(angle_attack);

        return cf;

    }

    double AerodynamicBody2D::incidentFlowPotential(Vector<double> x) const {
        return incident_flow(0) * x(0) + incident_flow(1) * x(1);
    }

}
