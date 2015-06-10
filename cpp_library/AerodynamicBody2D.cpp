/**
 * File:   AerodynamicBody2D.hpp
 * Author: kenobi
 *
 * Created on July 21, 2014, 9:29 AM
 */

#include"AerodynamicBody2D.hpp"
#include"Matrix.hpp"

#include <cmath>
#include <fstream>
#include <sstream>

namespace nde {

	AerodynamicBody2D::AerodynamicBody2D(
		   double chord,
		   const Vector<Panel2D>& panels)
	: _chord(chord), _panels(panels) {
		;
	}

	void AerodynamicBody2D::calcPotentialFlow(double angle_attack,
														PanelMethodType panel_method_type) {

		_angle_attack = angle_attack;
	  _incident_flow.resize(2);
	  _incident_flow(0) = cos(_angle_attack);
	  _incident_flow(1) = sin(_angle_attack);

		size_t num_panels = _panels.size();

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
				 for (size_t i = 0; i < num_panels; ++i) {
				     b(i) = _incident_flow * _panels(i).getControlPointIn()*(-1.0);
				     A(i, num_panels) = wake_panel.calcConstantDoubletPotencial
				             (_panels(i).getControlPointIn());
				     A(num_panels, i) = 0.0;
				     for (size_t j = 0; j < num_panels; ++j) {
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
				 for (size_t i = 0; i < num_panels; ++i)
				     doublets(i) = sol(i);
				 wake = sol(num_panels);

				/* calculate speed on the surface and local c_p*/
				_v_x.resize(num_panels - 1);
				_x.resize(num_panels - 1);
				_d_length_x.resize(num_panels - 1);
				_normal_x.resize(num_panels - 1);
				_cp_x.resize(num_panels - 1);
				for (size_t i = 0; i < num_panels - 1; ++i) {
			    	_x(i) = 0.5 * (_panels(i + 1).getMidPoint()
									 + _panels(i).getMidPoint()) / _chord;
					_d_length_x(i) = (_panels(i + 1).getMidPoint()
  								      - _panels(i).getMidPoint()).norm() / _chord;
			      _v_x(i) = std::abs((doublets(i + 1) - doublets(i))
										   / (_d_length_x(i) * _chord));
					_normal_x(i) = -0.5 * (_panels(i + 1).getNormal()
								 		      + _panels(i).getNormal());
					_cp_x(i) = 1 - _v_x(i) * _v_x(i);
				}

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
				 for (size_t i = 0; i < num_panels; ++i)
				     sources(i) = _incident_flow * _panels(i).getNormal();


				 /* calculate the doublets intensity */

				 // system matrix and independent vector
				 Matrix<double> A(num_panels + 1, num_panels + 1);
				 Vector<double> b(num_panels + 1);
				 for (size_t i = 0; i < num_panels; ++i) {
				     b(i) = 0.0;
				     A(i, num_panels) = wake_panel.calcConstantDoubletPotencial
				             (_panels(i).getControlPointIn());
				     A(num_panels, i) = 0.0;
				     for (size_t j = 0; j < num_panels; ++j) {
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
				 for (size_t i = 0; i < num_panels; ++i)
				     doublets(i) = sol(i);
				 wake = sol(num_panels);

				 /* calculate speed on the surface and local c_p*/
				_v_x.resize(num_panels - 1);
				_x.resize(num_panels - 1);
				_d_length_x.resize(num_panels - 1);
				_normal_x.resize(num_panels - 1);
				_cp_x.resize(num_panels - 1);
				for (size_t i = 0; i < num_panels - 1; ++i) {
 			      Vector<double> mt = 0.5 * (_panels(i + 1).getTangent()
													 + _panels(i).getTangent());
					_d_length_x(i) = (_panels(i + 1).getMidPoint()
  								      - _panels(i).getMidPoint()).norm() / _chord;
			      _v_x(i) = std::abs(_incident_flow * mt + (doublets(i + 1)
			 								 - doublets(i)) / (_d_length_x(i) * _chord));
			    	_x(i) = 0.5 * (_panels(i + 1).getMidPoint()
									 + _panels(i).getMidPoint()) / _chord;
					_normal_x(i) = -0.5 * (_panels(i + 1).getNormal()
								 		      + _panels(i).getNormal());
				   _cp_x(i) = 1 - _v_x(i) * _v_x(i);
				}

				break;

			}
			case NEUMANN_CONSTANT_SOURCES_AND_VORTEX:
			{

				 // system matrix and indep vector
				 Matrix<double> A(num_panels + 1, num_panels + 1);
				 Vector<double> b(num_panels + 1);

				 double sum_tangent_vortexes = 0.0;
				 for (size_t i = 0; i < num_panels; ++i) {

				     b(i) = _incident_flow * _panels(i).getNormal()*(-1.0);

				     double sum_normal_vortexes = 0.0;

				     for (size_t j = 0; j < num_panels; ++j) {
				         A(i, j) = _panels(j).calcConstantSourceSpeed
				                 (_panels(i).getControlPointOut())
									  * _panels(i).getNormal();
				         sum_normal_vortexes += _panels(j).calcConstantVortexSpeed
				                 (_panels(i).getControlPointOut())
									* _panels(i).getNormal();
				     }

				     A(i, num_panels) = sum_normal_vortexes;

				     A(num_panels, i) = _panels(i).calcConstantSourceSpeed
				             (_panels(0).getControlPointOut())
								* _panels(0).getTangent()
				             + _panels(i).calcConstantSourceSpeed
				             (_panels(num_panels - 1).getControlPointOut())
								* _panels(num_panels - 1).getTangent();

				     sum_tangent_vortexes += _panels(i).calcConstantVortexSpeed
				             (_panels(0).getControlPointOut())
								* _panels(0).getTangent()
				             + _panels(i).calcConstantVortexSpeed
				             (_panels(num_panels - 1).getControlPointOut())
								* _panels(num_panels - 1).getTangent();

				 }

				 A(num_panels, num_panels) = sum_tangent_vortexes;
				 b(num_panels) = -(_incident_flow * _panels(0).getTangent()
				         + _incident_flow * _panels(num_panels - 1).getTangent());

				 Vector<double> sources(num_panels);
				 double vortex;

				 Vector<double> sol = A.solve(b);

				 for (size_t i = 0; i < num_panels; ++i)
				     sources(i) = sol(i);
				 vortex = sol(num_panels);


				 /* calculate speed on the surface and local c_p*/
				_v_x.resize(num_panels);
				_x.resize(num_panels);
				_d_length_x.resize(num_panels);
				_normal_x.resize(num_panels);
				_cp_x.resize(num_panels);

				for (size_t i = 0; i < num_panels; ++i) {
					_x(i) = _panels(i).getControlPointOut() / _chord;
					Vector<double> p = _panels(i).getControlPointOut();
					_v_x(i) = _incident_flow * _panels(i).getTangent();
					for (size_t j = 0; j < num_panels; ++j)
					_v_x(i) += (_panels(j).calcConstantSourceSpeed(p) * sources(j)
						 + _panels(j).calcConstantVortexSpeed(p) * vortex)
	 					 * _panels(i).getTangent();
					_v_x(i) = std::abs(_v_x(i));
					_d_length_x(i) = _panels(i).getLength() / _chord;
					_normal_x(i) = -1.0 * _panels(i).getNormal();
					_cp_x(i) = 1 - _v_x(i) * _v_x(i);
				}

				break;

			}
			default:
				 //do nothing!
				 break;
		}

		/* calculate total aerodynamic force and moment */
		Vector<double> Fg(2);
		Fg.fill(0.0);
		for (size_t i = 0; i < _x.size(); ++i) {
			Fg = Fg + (_cp_x(i) * _d_length_x(i)) * _normal_x(i);
			_c_M0 = _c_M0 - (_cp_x(i) * _d_length_x(i)) * _x(i)(0) * _normal_x(i)(1);
		}

		_c_F.resize(2);
		_c_F(0) =  Fg(0) * cos(_angle_attack) + Fg(1) * sin(_angle_attack);
		_c_F(1) = -Fg(0) * sin(_angle_attack) + Fg(1) * cos(_angle_attack);

	}

	double AerodynamicBody2D::getAngleAttack() const {
		return _angle_attack;
	}

	Vector<double> AerodynamicBody2D::getForceCoeffs() const {
		return _c_F;
	}

	Vector<double> AerodynamicBody2D::getControlPointsDistances() const {

		size_t num_points = _x.size();
		Vector<double> d(num_points);
		d(0) = 0.0; // we take the first point as the origin to measure distances
		for (size_t i = 1; i < num_points; ++i)
			d(i) = d(i-1) + (_x(i)-_x(i-1)).norm(); // accumulated distance from
																 // the origin.
		return d;

	}

	Vector<double> AerodynamicBody2D::getVelocities() const {
		return _v_x;
	}

	void AerodynamicBody2D::writeResultsToFile
								(const std::string& file_name) const {

		std::ofstream out_file(file_name.c_str(), std::ofstream::out);

		out_file << "# d(x)" << "|" << "v(x)" << "|"
					<< "c_p(x)" << "|" << "dl(x)" << std::endl;

		Vector<double> d_x = getControlPointsDistances();
		for (size_t i = 0; i < _x.size(); ++i) {
			out_file << d_x(i) << " " << _v_x(i) << " "
						<< _cp_x(i) << " " << _d_length_x(i) << std::endl;
		}

	}


}

