/**
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
		DIRICHLET_CONSTANT_SOURCES_AND_DOUBLETS = 2,
		NEUMANN_CONSTANT_SOURCES_AND_VORTEX = 3
	};

	class AerodynamicBody2D {
	public:
		AerodynamicBody2D(
				 double chord,
				 const Vector<Panel2D>& panels);

		void calcPotentialFlow(double angle_attack,
									PanelMethodType panel_method_type);

		// property gets
		double getAngleAttack() const;
		Vector<double> getForceCoeffs() const;
		double getMomentCoeff() const;
		Vector<double> getControlPointsDistances() const;
		Vector<double> getVelocities() const;
		Vector<Vector<double> > getTangents() const;

		// write results
		void writeResultsToFile(const std::string& file_name) const;

	private:
		double _chord;
		double _angle_attack;
		Vector<Panel2D> _panels;
		Vector<double> _incident_flow;

		Vector<double> _v_x;
		Vector<Vector<double> > _x;
		Vector<double> _d_length_x;
		Vector<Vector<double> > _normal_x;
		Vector<Vector<double> > _tangent_x;
		Vector<double> _cp_x;

		Vector<double> _c_F;
		double _c_M0;

	};

}


#endif	/* AERODYNAMICBODY2D_HPP */

