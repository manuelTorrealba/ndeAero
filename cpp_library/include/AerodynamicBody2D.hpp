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
				 const Vector<Panel2D>& panels,
				 double angle_attack);

		void changeAngleAttack(double angle_attack);
		void calcPotentialFlow(PanelMethodType panel_method_type);
		Vector<double> getForceCoeffs() const;
		Vector<double> getControlPointsDistances(const Vector<double>& p0) const;
		Vector<double> getVelocities() const;


		private:
		double _chord;
		double _angle_attack;
		Vector<Panel2D> _panels;
		Vector<double> _incident_flow;

		Vector<double> _x;
		Vector<double> _v;
		Vector<double> _cp;
		Vector<double> _F;

	};

}


#endif	/* AERODYNAMICBODY2D_HPP */

