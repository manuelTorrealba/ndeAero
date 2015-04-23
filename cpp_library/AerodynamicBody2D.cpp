/**
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
		   double chord,
		   const Vector<Panel2D>& panels,
		   double angle_attack)
	: _chord(chord), _panels(panels), _angle_attack(angle_attack) {
	  _incident_flow.resize(2);
	  _incident_flow(0) = cos(_angle_attack);
	  _incident_flow(1) = sin(_angle_attack);
	}

	void AerodynamicBody2D::changeAngleAttack(double angle_attack) {
		_angle_attack = angle_attack;
	  _incident_flow.resize(2);
	  _incident_flow(0) = cos(_angle_attack);
	  _incident_flow(1) = sin(_angle_attack);
	}

    void AerodynamicBody2D::calcPotentialFlow
													(PanelMethodType panel_method_type) {

        int num_panels = _panels.size();

        switch (panel_method_type) {

            case DIRICHLET_CONSTANT_DOUBLETS:
            {

                Vector<double> far_wake_point = _panels(0).getStartPoint();
                far_wake_point(0) += 100.0;
                Panel2D wake_panel(_panels(0).getStartPoint(),
                        far_wake_point);

                /* calculate the doublets intensity */

                // system matrix and independent vector
                Matrix<double> A(num_panels + 1, num_panels + 1);
                Vector<double> b(num_panels + 1);
                for (int i = 0; i < num_panels; ++i) {
                    b(i) = _incident_flow * _panels(i).getControlPointIn()*(-1.0);
                    A(i, num_panels) = wake_panel.calcConstantDoubletPotencial
                            (_panels(i).getControlPointIn());
                    A(num_panels, i) = 0.0;
                    for (int j = 0; j < num_panels; ++j) {
                        A(i, j) = _panels(j).calcConstantDoubletPotencial
                                (_panels(i).getControlPointIn());
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
                _x.resize(num_panels - 1);
                _v.resize(num_panels - 1);
                _cp.resize(num_panels - 1);
                for (int i = 0; i < num_panels - 1; ++i) {
                    Vector<double> dl = _panels(i + 1).getMidPoint() 
												  - _panels(i).getMidPoint();
                    _x(i) = (_panels(i + 1).getMidPoint()(0)
								  + _panels(i).getMidPoint()(0)) * 0.5;
                    _v(i) = (doublets(i + 1) - doublets(i)) / dl.norm();
                    _cp(i) = 1 - _v(i) * _v(i);
                }

                /* calculate total aerodynamic force */
                Vector<double> Fg(2);
                Fg.fill(0.0);
                for (int i = 0; i < num_panels - 1; ++i) {
                    Vector<double> dl = _panels(i + 1).getMidPoint() 
												  - _panels(i).getMidPoint();
                    Vector<double> n = (_panels(i + 1).getNormal() 
												  + _panels(i).getNormal())* (-0.5);
                    Fg = Fg + n * (_cp(i) * dl.norm());
                }
                Fg = Fg / _chord;

                _F.resize(2);
                _F(0) = Fg(0) * cos(_angle_attack) + Fg(1) * sin(_angle_attack);
                _F(1) = -Fg(0) * sin(_angle_attack) + Fg(1) * cos(_angle_attack);

                break;

            }
            case DIRICHLET_CONSTANT_SOURCES_AND_DOUBLETS:
            {

                Vector<double> far_wake_point = _panels(0).getStartPoint();
                far_wake_point(0) += 100.0;
                Panel2D wake_panel(_panels(0).getStartPoint(),
                        far_wake_point);

                /* calculate the sources */
                Vector<double> sources(num_panels);
                for (int i = 0; i < num_panels; ++i)
                    sources(i) = _incident_flow * _panels(i).getNormal();


                /* calculate the doublets intensity */

                // system matrix and independent vector
                Matrix<double> A(num_panels + 1, num_panels + 1);
                Vector<double> b(num_panels + 1);
                for (int i = 0; i < num_panels; ++i) {
                    b(i) = 0.0;
                    A(i, num_panels) = wake_panel.calcConstantDoubletPotencial
                            (_panels(i).getControlPointIn());
                    A(num_panels, i) = 0.0;
                    for (int j = 0; j < num_panels; ++j) {
                        b(i) -= sources(j) * _panels(j).calcConstantSourcePotencial
                                (_panels(i).getControlPointIn());
                        A(i, j) = _panels(j).calcConstantDoubletPotencial
                                (_panels(i).getControlPointIn());
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
                _x.resize(num_panels - 1);
                _v.resize(num_panels - 1);
                _cp.resize(num_panels - 1);
                for (int i = 0; i < num_panels - 1; ++i) {
                    Vector<double> dl = _panels(i + 1).getMidPoint() 
												  - _panels(i).getMidPoint();
                    Vector<double> mt = (_panels(i + 1).getTangent() 
												   + _panels(i).getTangent()) * 0.5;
                    _x(i) = (_panels(i + 1).getMidPoint()(0) 
								  + _panels(i).getMidPoint()(0)) * 0.5;
                    _v(i) = _incident_flow * mt + (doublets(i + 1)
								 - doublets(i)) / dl.norm();
                    _cp(i) = 1 - _v(i) * _v(i);
                }

                /* calculate total aerodynamic force */
                Vector<double> Fg(2);
                Fg.fill(0.0);
                for (int i = 0; i < num_panels - 1; ++i) {
                    Vector<double> dl = _panels(i + 1).getMidPoint() 
												  - _panels(i).getMidPoint();
                    Vector<double> n = (_panels(i + 1).getNormal() 
												 + _panels(i).getNormal())* (-0.5);
                    Fg = Fg + n * (_cp(i) * dl.norm());
                }
                Fg = Fg / _chord;

                _F.resize(2);
                _F(0) = Fg(0) * cos(_angle_attack) + Fg(1) * sin(_angle_attack);
                _F(1) = -Fg(0) * sin(_angle_attack) + Fg(1) * cos(_angle_attack);

                break;

            }
            case NEUMANN_CONSTANT_SOURCES_AND_VORTEX:
            {

                // system matrix and indep vector
                Matrix<double> A(num_panels + 1, num_panels + 1);
                Vector<double> b(num_panels + 1);

                double sum_tangent_vortexes = 0.0;
                for (int i = 0; i < num_panels; ++i) {

                    b(i) = _incident_flow * _panels(i).getNormal()*(-1.0);

                    double sum_normal_vortexes = 0.0;

                    for (int j = 0; j < num_panels; ++j) {
                        A(i, j) = _panels(j).calcConstantSourceSpeed
                                (_panels(i).getControlPointOut()) * _panels(i).getNormal();
                        sum_normal_vortexes += _panels(j).calcConstantVortexSpeed
                                (_panels(i).getControlPointOut()) * _panels(i).getNormal();
                    }

                    A(i, num_panels) = sum_normal_vortexes;

                    A(num_panels, i) = _panels(i).calcConstantSourceSpeed
                            (_panels(0).getControlPointOut()) * _panels(0).getTangent()
                            + _panels(i).calcConstantSourceSpeed
                            (_panels(num_panels - 1).getControlPointOut()) * _panels(num_panels - 1).getTangent();

                    sum_tangent_vortexes += _panels(i).calcConstantVortexSpeed
                            (_panels(0).getControlPointOut()) * _panels(0).getTangent()
                            + _panels(i).calcConstantVortexSpeed
                            (_panels(num_panels - 1).getControlPointOut()) * _panels(num_panels - 1).getTangent();

                }

                A(num_panels, num_panels) = sum_tangent_vortexes;
                b(num_panels) = -(_incident_flow * _panels(0).getTangent()
                        + _incident_flow * _panels(num_panels - 1).getTangent());

                Vector<double> sources(num_panels);
                double vortex;

                Vector<double> sol = A.solve(b);

                for (int i = 0; i < num_panels; ++i)
                    sources(i) = sol(i);
                vortex = sol(num_panels);


                /* calculate speed on the surface and local c_p*/
                _x.resize(num_panels - 1);
                _v.resize(num_panels - 1);
                _cp.resize(num_panels - 1);
                Vector<double> Fg(2);
                Fg.fill(0.0);
                for (int i = 0; i < num_panels - 1; ++i) {
                    _x(i) = _panels(i).getControlPointOut()(0);
                    Vector<double> p = _panels(i).getControlPointOut();
                    _v(i) = _incident_flow * _panels(i).getTangent();
                    for (int j = 0; j < num_panels; ++j)
                        _v(i) += (_panels(j).calcConstantSourceSpeed(p) * sources(j)
                            + _panels(j).calcConstantVortexSpeed(p) * vortex)
                        * _panels(i).getTangent();
                    _cp(i) = 1 - _v(i) * _v(i);
                    Fg = Fg - _panels(i).getNormal() * _cp(i) * _panels(i).getLength();
                }
                Fg = Fg / _chord;

                _F.resize(2);
                _F(0) = Fg(0) * cos(_angle_attack) + Fg(1) * sin(_angle_attack);
                _F(1) = -Fg(0) * sin(_angle_attack) + Fg(1) * cos(_angle_attack);

                break;

            }

            default:
                //do nothing!
                break;
        }

    }

    Vector<double> AerodynamicBody2D::getForceCoeffs() const {
        return _F;
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
