/* 
 * File:   Panel.cpp
 * Author: kenobi
 *
 * Created on July 21, 2014, 9:34 AM
 */

#include"Panel.hpp"
#include"MathUtils.hpp"

#include<cmath>

using namespace std;

namespace nde {

/*********************************************************************************************/
/*
 * 2D PANEL CLASS
 */
/*********************************************************************************************/

Panel2D::Panel2D(const Vector<double>& start_point_in,
		const Vector<double>& end_point_in) :
		start_point(start_point_in), end_point(end_point_in) {

	setPanel2D();

}

void Panel2D::setPoints(const Vector<double>& start_point_in,
		const Vector<double>& end_point_in) {

	start_point = start_point_in;
	end_point = end_point_in;
	setPanel2D();

}

void Panel2D::setPanel2D() {

	dl = end_point - start_point;
	length = dl.norm();
	tangent = dl / length;
	mid_point = (start_point + end_point) * 0.5;
	normal.resize(2);
	normal(0) = -tangent(1);
	normal(1) = tangent(0);

	double theta = atan2(tangent(1), tangent(0));
	B.resize(2, 2);
	B(0, 0) = cos(theta);
	B(0, 1) = sin(theta);
	B(1, 0) = -sin(theta);
	B(1, 1) = cos(theta);

	start_point_local = B * start_point;
	end_point_local = B * end_point;

	control_point_in = mid_point - normal * (length * 1.0e-4);
	control_point_out = mid_point + normal * (length * 1.0e-4);

}

Vector<double> Panel2D::getStartPoint() const {
	return start_point;
}

Vector<double> Panel2D::getEndPoint() const {
	return end_point;
}

double Panel2D::getDx() const {
	return dl(0);
}

double Panel2D::getDz() const {
	return dl(1);
}

Vector<double> Panel2D::getDl() const {
	return dl;
}

double Panel2D::getLength() const {
	return length;
}

Vector<double> Panel2D::getMidPoint() const {
	return mid_point;
}

Vector<double> Panel2D::getControlPointIn() const {
	return control_point_in;
}

Vector<double> Panel2D::getControlPointOut() const {
	return control_point_out;
}

Vector<double> Panel2D::getNormal() const {
	return normal;
}

Vector<double> Panel2D::getTangent() const {
	return tangent;
}

double Panel2D::calcConstantDoubletPotencial(const Vector<double>& x) const {

	Vector<double> xl = B * x;

	double theta1 = atan2(xl(1) - start_point_local(1),
			xl(0) - start_point_local(0));
	double theta2 = atan2(xl(1) - end_point_local(1),
			xl(0) - end_point_local(0));

	return -(theta2 - theta1) / (2.0 * M_PI);

}

double Panel2D::calcConstantSourcePotencial(const Vector<double>& x) const {

	double r1 = (x - start_point).norm();
	double r2 = (x - end_point).norm();
	Vector<double> xl = B * x;
	double theta1 = atan2(xl(1) - start_point_local(1),
			xl(0) - start_point_local(0));
	double theta2 = atan2(xl(1) - end_point_local(1),
			xl(0) - end_point_local(0));

	double tol = 1.0e-6 * length;

	double t1;
	if (r1 > tol)
		t1 = 2.0 * r1 * cos(theta1) * log(r1);
	else
		t1 = 0.0;

	double t2;
	if (r2 > tol)
		t2 = -2.0 * r2 * cos(theta2) * log(r2);
	else
		t2 = 0.0;

	double t3 = 2.0 * r1 * sin(theta1) * (theta2 - theta1);

	return (t1 + t2 + t3) / (4.0 * M_PI);

}

Vector<double> Panel2D::calcConstantDoubletSpeed(
		const Vector<double>& x) const {

	Vector<double> u(2);
	Vector<double> r1 = x - start_point;
	Vector<double> r2 = x - end_point;
	double r1sqr = r1.norm() * r1.norm();
	double r2sqr = r2.norm() * r2.norm();

	u(0) = 1.0 / (2.0 * M_PI) * (r2(1) / r2sqr - r1(1) / r1sqr);
	u(1) = -1.0 / (2.0 * M_PI) * (r2(0) / r2sqr - r1(0) / r1sqr);

	return u;

}

Vector<double> Panel2D::calcConstantSourceSpeed(const Vector<double>& x) const {

	Vector<double> xl = B * x;
	double theta1 = atan2(xl(1) - start_point_local(1),
			xl(0) - start_point_local(0));
	double theta2 = atan2(xl(1) - end_point_local(1),
			xl(0) - end_point_local(0));
	Vector<double> r1 = x - start_point;
	Vector<double> r2 = x - end_point;

	Vector<double> u(2);
	u(0) = 1.0 / (2.0 * M_PI) * log(r1.norm() / r2.norm());
	u(1) = 1.0 / (2.0 * M_PI) * (theta2 - theta1);

	return B.solve(u); //return in global coordinates

}

Vector<double> Panel2D::calcConstantVortexSpeed(const Vector<double>& x) const {

	Vector<double> xl = B * x;
	double theta1 = atan2(xl(1) - start_point_local(1),
			xl(0) - start_point_local(0));
	double theta2 = atan2(xl(1) - end_point_local(1),
			xl(0) - end_point_local(0));
	Vector<double> r1 = x - start_point;
	Vector<double> r2 = x - end_point;

	Vector<double> u(2);
	u(0) = 1.0 / (2.0 * M_PI) * (theta2 - theta1);
	u(1) = 1.0 / (2.0 * M_PI) * log(r2.norm() / r1.norm());

	return B.solve(u); //return in global coordinates

}

/*********************************************************************************************/
/*
 * 3D PANEL CLASS
 */
/*********************************************************************************************/

Panel3D::Panel3D(const Vector<double>& x1_in, const Vector<double>& x2_in,
		const Vector<double>& x3_in, const Vector<double>& x4_in,
		const Vector<double>& xc_in) :
		x1(x1_in), x2(x2_in), x3(x3_in), x4(x4_in), xc(xc_in) {
	;
}

void Panel3D::setPoints(const Vector<double>& x1_in,
		const Vector<double>& x2_in, const Vector<double>& x3_in,
		const Vector<double>& x4_in, const Vector<double>& xc_in) {

	x1 = x1_in;
	x2 = x2_in;
	x3 = x3_in;
	x4 = x4_in;
	xc = xc_in;

}

void Panel3D::setSourceIntensities(double sc12c_in, double sc23c_in,
		double sc34c_in, double sc14c_in) {
	sc12c = sc12c_in;
	sc23c = sc23c_in;
	sc34c = sc34c_in;
	sc14c = sc14c_in;
}

void Panel3D::setDoubletIntensities(double dt12c_in, double dt23c_in,
		double dt34c_in, double dt14c_in) {
	dt12c = dt12c_in;
	dt23c = dt23c_in;
	dt34c = dt34c_in;
	dt14c = dt14c_in;
}

double Panel3D::getArea() const {
	double area = 0.0;
	for (int i = 0; i < 4; ++i)
		area += sub_panels(i).getArea();
	return area;
}

void Panel3D::setPanel3D() {

	/* General geometry of the panels ->
	 * Square panel composed of 4 triangular panels
	 *         +(1)+++++++++++(4)+
	 *		   +   \                        /     +
	 *         +      \                  /        +
	 *         +         \            /           +
	 *         +            +(c)               +
	 *         +          /           \           +
	 *         +       /                 \        +
	 *         +    /                       \     +
	 *         +(2)+++++++++++(3)+
	 */

	sub_panels.resize(4);
	// triangle 12c
	sub_panels(0).setPoints(x1, x2, xc);
	// triangle 23c
	sub_panels(1).setPoints(x2, x3, xc);
	// triangle 34c
	sub_panels(2).setPoints(x3, x4, xc);
	// triangle 41c
	sub_panels(3).setPoints(x4, x1, xc);

}

// access a subpanel
const SubPanel3DTriangle& Panel3D::getSubPanel(int index) const {
	return sub_panels(index);
}

double Panel3D::getVelocity() const {

	// Solve for the polynomial coefficients for the doublets intensities.
	// Assume doublet is applied at the centroid of each triangular panel.
	// Functional form:
	// y = coeff(0) + coeff(1) x + coeff(2) y + coeff(3) x y
	// Fit coeff(i) to the four triangular sub-panels

	// center of masses
	Vector<double> cm_12c = sub_panels(0).getCenterOfMass();
	Vector<double> cm_23c = sub_panels(1).getCenterOfMass();
	Vector<double> cm_34c = sub_panels(2).getCenterOfMass();
	Vector<double> cm_41c = sub_panels(3).getCenterOfMass();

	// set up equations
	Matrix<double> A(4, 4);
	Vector<double> b(4);

	A(0, 0) = 1.0;
	A(0, 1) = cm_12c(0);
	A(0, 2) = cm_12c(1);
	A(0, 3) = cm_12c(0) * cm_12c(1);
	b(0) = dt12c;

	A(1, 0) = 1.0;
	A(1, 1) = cm_23c(0);
	A(1, 2) = cm_23c(1);
	A(1, 3) = cm_23c(0) * cm_23c(1);
	b(1) = dt23c;

	A(2, 0) = 1.0;
	A(2, 1) = cm_34c(0);
	A(2, 2) = cm_34c(1);
	A(2, 3) = cm_34c(0) * cm_34c(1);
	b(2) = dt34c;

	A(3, 0) = 1.0;
	A(3, 1) = cm_41c(0);
	A(3, 2) = cm_41c(1);
	A(3, 3) = cm_41c(0) * cm_41c(1);
	b(3) = dt14c;

	Vector<double> coeff = A.solve(b);

	// Velocity is given by derivative of doublet intensity along x and y directions
	double v_x = coeff(1) + coeff(3) * xc(1);
	double v_y = coeff(2) + coeff(3) * xc(0);

	return sqrt(v_x * v_x + v_y * v_y);

}

Vector<double> Panel3D::getForce() const {

	double v = getVelocity();
	double p = 1.0 - 0.5 * v * v;

	Vector<double> f(3);
	f.fill(0.0);
	for (int i = 0; i < 4; ++i)
		f = f + sub_panels(i).getNormal() * sub_panels(i).getArea();

	return f * p;

}

Vector<double> Panel3D::getMoment(const Vector<double>& x) const {

	double v = getVelocity();
	double p = 1.0 - 0.5 * v * v;

	Vector<double> moment(3);
	moment.fill(0.0);
	for (int i = 0; i < 4; ++i) {
		moment = moment
				+ calcCrossProduct(sub_panels(i).getCenterOfMass() - x,
						sub_panels(i).getNormal())
						* (p * sub_panels(i).getArea());
	}

	return moment;

}

SubPanel3DTriangle::SubPanel3DTriangle(const Vector<double>& x1_in,
		const Vector<double>& x2_in, const Vector<double>& x3_in) :
		x1(x1_in), x2(x2_in), x3(x3_in) {

	setSubPanel3DTriangle();

}

void SubPanel3DTriangle::setPoints(const Vector<double>& x1_in,
		const Vector<double>& x2_in, const Vector<double>& x3_in) {

	x1 = x1_in;
	x2 = x2_in;
	x3 = x3_in;

	setSubPanel3DTriangle();

}

void SubPanel3DTriangle::setSubPanel3DTriangle() {

	Vector<double> x12 = x2 - x1;
	Vector<double> x13 = x3 - x1;
	Vector<double> xx123 = calcCrossProduct(x12, x13);

	area = 0.5 * xx123.norm();
	normal = xx123 / (2.0 * area);
	center_of_mass = (x1 + x2 + x3) / 3.0;

	//TODO: define matrix B -> from global to local coordinates!
	x1l = B * (x1 - center_of_mass);
	x2l = B * (x2 - center_of_mass);
	x3l = B * (x3 - center_of_mass);

}

double SubPanel3DTriangle::getArea() const {
	return area;
}

Vector<double> SubPanel3DTriangle::getNormal() const {
	return normal;
}

Vector<double> SubPanel3DTriangle::getCenterOfMass() const {
	return center_of_mass;
}

double SubPanel3DTriangle::calcConstantSourcePotential(
		const Vector<double>& x) const {

	Vector<double> xl = B * (x - center_of_mass);

	double z = xl(2);

	double r1 = (xl - x1l).norm();
	double r2 = (xl - x2l).norm();
	double r3 = (xl - x3l).norm();

	double d12 = (x2l - x1l).norm();
	double d23 = (x3l - x2l).norm();
	double d31 = (x1l - x3l).norm();

	double m12 = (x2l(1) - x1l(1)) / (x2l(0) - x1l(0));
	double m23 = (x3l(1) - x2l(1)) / (x3l(0) - x2l(0));
	double m31 = (x1l(1) - x3l(1)) / (x1l(0) - x3l(0));

	double e1 = pow(xl(0) - x1l(0), 2) + pow(xl(2) - x1l(2), 2);
	double e2 = pow(xl(0) - x2l(0), 2) + pow(xl(2) - x2l(2), 2);
	double e3 = pow(xl(0) - x3l(0), 2) + pow(xl(2) - x3l(2), 2);

	double h1 = (xl(0) - x1l(0)) * (xl(1) - x1l(1));
	double h2 = (xl(0) - x2l(0)) * (xl(1) - x2l(1));
	double h3 = (xl(0) - x3l(0)) * (xl(1) - x3l(1));

	double mln12 = ((xl(0) - x1l(0)) * (x2l(1) - x1l(1))
			- (xl(1) - x1l(1)) * (x2l(0) - x1l(0))) / d12;
	double mln23 = ((xl(0) - x2l(0)) * (x3l(1) - x2l(1))
			- (xl(1) - x2l(1)) * (x3l(0) - x2l(0))) / d23;
	double mln31 = ((xl(0) - x3l(0)) * (x1l(1) - x3l(1))
			- (xl(1) - x3l(1)) * (x1l(0) - x3l(0))) / d31;

	double at12 = atan2(m12 * e1 - h1, z * r1) - atan2(m12 * e2 - h2, z * r2);
	double at23 = atan2(m23 * e2 - h2, z * r2) - atan2(m23 * e3 - h3, z * r3);
	double at31 = atan2(m31 * e3 - h3, z * r3) - atan2(m31 * e1 - h1, z * r1);

	double ln12 = log((r1 + r2 + d12) / (r1 + r2 - d12));
	double ln23 = log((r2 + r3 + d23) / (r2 + r3 - d23));
	double ln31 = log((r3 + r1 + d31) / (r3 + r1 - d31));

	return -(mln12 * ln12 + mln23 * ln23 + mln31 * ln31
			+ abs(z) * (at12 + at23 + at31)) / (4.0 * M_PI);

}

double SubPanel3DTriangle::calcConstantDoubletPotential(
		const Vector<double>& x) const {

	Vector<double> xl = B * (x - center_of_mass);

	double z = xl(2);

	double r1 = (xl - x1l).norm();
	double r2 = (xl - x2l).norm();
	double r3 = (xl - x3l).norm();

	double d12 = (x2l - x1l).norm();
	double d23 = (x3l - x2l).norm();
	double d31 = (x1l - x3l).norm();

	double m12 = (x2l(1) - x1l(1)) / (x2l(0) - x1l(0));
	double m23 = (x3l(1) - x2l(1)) / (x3l(0) - x2l(0));
	double m31 = (x1l(1) - x3l(1)) / (x1l(0) - x3l(0));

	double e1 = pow(xl(0) - x1l(0), 2) + pow(xl(2) - x1l(2), 2);
	double e2 = pow(xl(0) - x2l(0), 2) + pow(xl(2) - x2l(2), 2);
	double e3 = pow(xl(0) - x3l(0), 2) + pow(xl(2) - x3l(2), 2);

	double h1 = (xl(0) - x1l(0)) * (xl(1) - x1l(1));
	double h2 = (xl(0) - x2l(0)) * (xl(1) - x2l(1));
	double h3 = (xl(0) - x3l(0)) * (xl(1) - x3l(1));

	double at12 = atan2(m12 * e1 - h1, z * r1) - atan2(m12 * e2 - h2, z * r2);
	double at23 = atan2(m23 * e2 - h2, z * r2) - atan2(m23 * e3 - h3, z * r3);
	double at31 = atan2(m31 * e3 - h3, z * r3) - atan2(m31 * e1 - h1, z * r1);

	return (at12 + at23 + at31) / (4.0 * M_PI);

}

} /*end namespace nde*/
