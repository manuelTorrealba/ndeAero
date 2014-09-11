/*
 * File: MathUtils.cpp
 * Author: kenobi
 *
 * Created on September 9, 2014, 15:55
 */

#include"MathUtils.hpp"
#include"Vector.hpp"

namespace nde {

Vector<double> calcCrossProduct(const Vector<double>& v1,
		const Vector<double>& v2) {

	Vector<double> w(3);

	w(0) = v1(1) * v2(2) - v2(1) * v1(2);
	w(1) = -(v1(0) * v2(2) - v2(0) * v1(2));
	w(2) = v1(0) * v2(1) - v2(0) * v1(1);

	return w;

}

} /*end of namespace nde */
