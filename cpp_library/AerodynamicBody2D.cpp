/* 
 * File:   AerodynamicBody2D.hpp
 * Author: kenobi
 *
 * Created on July 21, 2014, 9:29 AM
 */

#include"AerodynamicBody2D.hpp"
#include"Matrix.hpp"

#include<cmath>

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

    void AerodynamicBody2D::calcPotentialFlow(PanelMethodType panel_method_type) {

        int num_panels = panels.size();

        switch (panel_method_type) {

            case DIRICHLET_CONSTANT_DOUBLETS:
            {

                Vector<double> far_wake_point = panels(0).getStartPoint();
                far_wake_point(0) += 100.0;
                Panel2D wake_panel(panels(0).getStartPoint(),
                        far_wake_point);

                /* calculate the doublets intensity */

                // system matrix and independent vector
                Matrix<double> A(num_panels + 1, num_panels + 1);
                Vector<double> b(num_panels + 1);
                for (int i = 0; i < num_panels; ++i) {
                    b(i) = incident_flow * panels(i).getControlPointIn()*(-1.0);
                    A(i, num_panels) = wake_panel.calcConstantDoubletPotencial
                            (panels(i).getControlPointIn());
                    A(num_panels, i) = 0.0;
                    for (int j = 0; j < num_panels; ++j) {
                        A(i, j) = panels(j).calcConstantDoubletPotencial
                                (panels(i).getControlPointIn());
                    }
                }
                A(num_panels, 0) = 1.0;
                A(num_panels, num_panels - 1) = -1.0;
                A(num_panels, num_panels) = 1.0;
                b(num_panels) = 0.0;

                Vector<double> sol = A.solve(b);

                Vector<double> doublets(num_panels);
                double wake;
                for (int i = 0; i < num_panels; ++i)
                    doublets(i) = sol(i);
                wake = sol(num_panels);

                /* calculate speed on the surface and local c_p*/
                x.resize(num_panels - 1);
                v.resize(num_panels - 1);
                cp.resize(num_panels - 1);
                for (int i = 0; i < num_panels - 1; ++i) {
                    Vector<double> dl = panels(i + 1).getMidPoint() - panels(i).getMidPoint();
                    x(i) = (panels(i + 1).getMidPoint()(0) + panels(i).getMidPoint()(0)) * 0.5;
                    v(i) = (doublets(i + 1) - doublets(i)) / dl.norm();
                    cp(i) = 1 - v(i) * v(i);
                }

                /* calculate total aerodynamic force */
                Vector<double> Fg(2);
                Fg.fill(0.0);
                for (int i = 0; i < num_panels - 1; ++i) {
                    Vector<double> dl = panels(i + 1).getMidPoint() - panels(i).getMidPoint();
                    Vector<double> n = (panels(i + 1).getNormal() + panels(i).getNormal())* (-0.5);
                    Fg = Fg + n * (cp(i) * dl.norm());
                }
                Fg = Fg / chord;

                F.resize(2);
                F(0) = Fg(0) * cos(angle_attack) + Fg(1) * sin(angle_attack);
                F(1) = -Fg(0) * sin(angle_attack) + Fg(1) * cos(angle_attack);

                break;

            }
            case DIRICHLET_CONSTANT_SOURCES_AND_DOUBLETS:
            {

                Vector<double> far_wake_point = panels(0).getStartPoint();
                far_wake_point(0) += 100.0;
                Panel2D wake_panel(panels(0).getStartPoint(),
                        far_wake_point);

                /* calculate the sources */
                Vector<double> sources(num_panels);
                for (int i = 0; i < num_panels; ++i)
                    sources(i) = incident_flow * panels(i).getNormal();


                /* calculate the doublets intensity */

                // system matrix and independent vector
                Matrix<double> A(num_panels + 1, num_panels + 1);
                Vector<double> b(num_panels + 1);
                for (int i = 0; i < num_panels; ++i) {
                    b(i) = 0.0;
                    A(i, num_panels) = wake_panel.calcConstantDoubletPotencial
                            (panels(i).getControlPointIn());
                    A(num_panels, i) = 0.0;
                    for (int j = 0; j < num_panels; ++j) {
                        b(i) -= sources(j) * panels(j).calcConstantSourcePotencial
                                (panels(i).getControlPointIn());
                        A(i, j) = panels(j).calcConstantDoubletPotencial
                                (panels(i).getControlPointIn());
                    }
                }
                A(num_panels, 0) = 1.0;
                A(num_panels, num_panels - 1) = -1.0;
                A(num_panels, num_panels) = 1.0;
                b(num_panels) = 0.0;

                Vector<double> sol = A.solve(b);

                Vector<double> doublets(num_panels);
                double wake;
                for (int i = 0; i < num_panels; ++i)
                    doublets(i) = sol(i);
                wake = sol(num_panels);

                /* calculate speed on the surface and local c_p*/
                x.resize(num_panels - 1);
                v.resize(num_panels - 1);
                cp.resize(num_panels - 1);
                for (int i = 0; i < num_panels - 1; ++i) {
                    Vector<double> dl = panels(i + 1).getMidPoint() - panels(i).getMidPoint();
                    Vector<double> mt = (panels(i + 1).getTangent() + panels(i).getTangent())*0.5;
                    x(i) = (panels(i + 1).getMidPoint()(0) + panels(i).getMidPoint()(0)) * 0.5;
                    v(i) = incident_flow * mt + (doublets(i + 1) - doublets(i)) / dl.norm();
                    cp(i) = 1 - v(i) * v(i);
                }

                /* calculate total aerodynamic force */
                Vector<double> Fg(2);
                Fg.fill(0.0);
                for (int i = 0; i < num_panels - 1; ++i) {
                    Vector<double> dl = panels(i + 1).getMidPoint() - panels(i).getMidPoint();
                    Vector<double> n = (panels(i + 1).getNormal() + panels(i).getNormal())* (-0.5);
                    Fg = Fg + n * (cp(i) * dl.norm());
                }
                Fg = Fg / chord;

                F.resize(2);
                F(0) = Fg(0) * cos(angle_attack) + Fg(1) * sin(angle_attack);
                F(1) = -Fg(0) * sin(angle_attack) + Fg(1) * cos(angle_attack);

                break;
            }
            case NEUMANN_CONSTANT_SOURCES_AND_VORTEX:
            {

                // system matrix and indep vector
                Matrix<double> A(num_panels + 1, num_panels + 1);
                Vector<double> b(num_panels + 1);

                double sum_tangent_vortexes = 0.0;
                for (int i = 0; i < num_panels; ++i) {

                    b(i) = incident_flow * panels(i).getNormal()*(-1.0);

                    double sum_normal_vortexes = 0.0;

                    for (int j = 0; j < num_panels; ++j) {
                        A(i, j) = panels(j).calcConstantSourceSpeed
                                (panels(i).getControlPointOut()) * panels(i).getNormal();
                        sum_normal_vortexes += panels(j).calcConstantVortexSpeed
                                (panels(i).getControlPointOut()) * panels(i).getNormal();
                    }

                    A(i, num_panels) = sum_normal_vortexes;

                    A(num_panels, i) = panels(i).calcConstantSourceSpeed
                            (panels(0).getControlPointOut()) * panels(0).getTangent()
                            + panels(i).calcConstantSourceSpeed
                            (panels(num_panels - 1).getControlPointOut()) * panels(num_panels - 1).getTangent();

                    sum_tangent_vortexes += panels(i).calcConstantVortexSpeed
                            (panels(0).getControlPointOut()) * panels(0).getTangent()
                            + panels(i).calcConstantVortexSpeed
                            (panels(num_panels - 1).getControlPointOut()) * panels(num_panels - 1).getTangent();

                }

                A(num_panels, num_panels) = sum_tangent_vortexes;
                b(num_panels) = -(incident_flow * panels(0).getTangent()
                        + incident_flow * panels(num_panels - 1).getTangent());

                Vector<double> sources(num_panels);
                double vortex;

                Vector<double> sol = A.solve(b);

                for (int i = 0; i < num_panels; ++i)
                    sources(i) = sol(i);
                vortex = sol(num_panels);


                /* calculate speed on the surface and local c_p*/
                x.resize(num_panels - 1);
                v.resize(num_panels - 1);
                cp.resize(num_panels - 1);
                Vector<double> Fg(2);
                Fg.fill(0.0);
                for (int i = 0; i < num_panels - 1; ++i) {
                    x(i) = panels(i).getControlPointOut()(0);
                    Vector<double> p = panels(i).getControlPointOut();
                    v(i) = incident_flow * panels(i).getTangent();
                    for (int j = 0; j < num_panels; ++j)
                        v(i) += (panels(j).calcConstantSourceSpeed(p) * sources(j)
                            + panels(j).calcConstantVortexSpeed(p) * vortex)
                        * panels(i).getTangent();
                    cp(i) = 1 - v(i) * v(i);
                    Fg = Fg - panels(i).getNormal() * cp(i) * panels(i).getLength();
                }
                Fg = Fg / chord;

                F.resize(2);
                F(0) = Fg(0) * cos(angle_attack) + Fg(1) * sin(angle_attack);
                F(1) = -Fg(0) * sin(angle_attack) + Fg(1) * cos(angle_attack);

                break;

            }

            default:
                //do nothing!
                break;
        }

    }

    Vector<double> AerodynamicBody2D::getForceCoeffs() const {
        return F;
    }

    //    Vector<double> AerodynamicBody2D::calcForceCoefficients() const {
    //
    //        Vector<double> cf_global(2);
    //        cf_global.fill(0.0);
    //
    //        for (int i = 0; i < panels.size() - 1; ++i) {
    //            Vector<double> v = getSpeed(panels(i).getControlPointOut());
    //            cf_global = cf_global - panels(i).getNormal()
    //                    * ((1.0 - pow(v.norm(), 2))
    //                    * panels(i).getLength() / chord);
    //        }
    //
    //        Vector<double> cf(2);
    //        cf(0) = cf_global(0) * cos(angle_attack) + cf_global(1) * sin(angle_attack);
    //        cf(1) = -cf_global(0) * sin(angle_attack) + cf_global(1) * cos(angle_attack);
    //
    //        return cf;
    //
    //    }


}