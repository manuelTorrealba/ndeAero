/* 
 * File:   Panel.hpp
 * Author: kenobi
 *
 * Created on July 21, 2014, 9:34 AM
 */

#ifndef PANEL_HPP
#define	PANEL_HPP

#include"Vector.hpp"
#include"Matrix.hpp"

namespace nde {

class Panel2D {
public:
	Panel2D(const Vector<double>& start_point_in,
			const Vector<double>& end_point_in);

	void setPoints(const Vector<double>& start_point_in,
			const Vector<double>& end_point_in);

	Vector<double> getStartPoint() const;
	Vector<double> getEndPoint() const;
	double getDx() const;
	double getDz() const;
	Vector<double> getDl() const;
	double getLength() const;
	Vector<double> getMidPoint() const;
	Vector<double> getControlPointIn() const;
	Vector<double> getControlPointOut() const;
	Vector<double> getNormal() const;
	Vector<double> getTangent() const;

	double calcConstantDoubletPotencial(const Vector<double>& x) const;
	double calcConstantSourcePotencial(const Vector<double>& x) const;

	Vector<double> calcConstantSourceSpeed(const Vector<double>& x) const;
	Vector<double> calcConstantDoubletSpeed(const Vector<double>& x) const;
	Vector<double> calcConstantVortexSpeed(const Vector<double>& x) const;

private:
	Vector<double> start_point;
	Vector<double> end_point;
	Vector<double> start_point_local;
	Vector<double> end_point_local;
	Vector<double> dl;
	double length;
	Vector<double> mid_point;
	Vector<double> control_point_in;
	Vector<double> control_point_out;
	Vector<double> normal;
	Vector<double> tangent;
	Matrix<double> B; //matrix to change from global to local coordinates
	void setPanel2D();

};

class SubPanel3DTriangle {
public:
	SubPanel3DTriangle(const Vector<double>& x1_in, const Vector<double>& x2_in,
			const Vector<double>& x3_in);

	void setPoints(const Vector<double>& x1_in, const Vector<double>& x2_in,
			const Vector<double>& x3_in);

	double getArea() const;
	Vector<double> getNormal() const;
	Vector<double> getCenterOfMass() const;

	double calcConstantSourcePotential(const Vector<double>& x) const;
	double calcConstantDoubletPotential(const Vector<double>& x) const;

private:
	Vector<double> x1, x2, x3;
	Vector<double> x1l, x2l, x3l;
	Matrix<double> B; //change from global to local coordinates matrix.
	Vector<double> normal;
	Vector<double> center_of_mass;
	double area;

	void setSubPanel3DTriangle();

};

class Panel3D {
public:

	// panel defined by four corner points.
	Panel3D(const Vector<double>& x1_in, const Vector<double>& x2_in,
			const Vector<double>& x3_in, const Vector<double>& x4_in,
			const Vector<double>& xc_in);

	void setPoints(const Vector<double>& x1_in, const Vector<double>& x2_in,
			const Vector<double>& x3_in, const Vector<double>& x4_in,
			const Vector<double>& xc_in);

	double getArea() const;

	void setSourceIntensities(double sc12c_in, double sc23c_in, double sc34c_in,
			double sc14c_in);
	void setDoubletIntensities(double dt12c_in, double dt23c_in,
			double dt34c_in, double dt14c_in);

	// access a subpanel
	const SubPanel3DTriangle& getSubPanel(int index) const;

	double getVelocity() const;
	Vector<double> getForce() const;
	Vector<double> getMoment(const Vector<double>& xm) const;

	double getSurface() const;
	Vector<double> getMidPoint() const;
	Vector<double> getControlPointIn() const;
	Vector<double> getControlPointOut() const;
	Vector<double> getNormal() const;

private:
	// constructor variables

	// panel corners
	Vector<double> x1;
	Vector<double> x2;
	Vector<double> x3;
	Vector<double> x4;
	Vector<double> xc; // control point

	// source intensities
	double sc12c;
	double sc23c;
	double sc34c;
	double sc14c;
	// doublet intensities
	double dt12c;
	double dt23c;
	double dt34c;
	double dt14c;

	/* derived variables */
	Vector<SubPanel3DTriangle> sub_panels;

	/* private methods */
	void setPanel3D();

};

}/*end of namespace nde */

#endif	/* PANEL_HPP */

