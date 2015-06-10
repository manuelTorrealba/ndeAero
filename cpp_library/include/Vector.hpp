/**
 * File:   PotentialFlowElements.hpp
 * Author: kenobi
 *
 * Created on July 17, 2014, 11:09 AM
 */

#ifndef VECTOR_HPP
#define VECTOR_HPP

#include<cstdlib>
#include<iostream>
#include<cmath>

namespace nde {

/* the template class declaration */
template<class T>
class Vector;


/* multiplication operators */
template <class U>
U operator*(const Vector<U>& lhs, const Vector<U>& rhs) {
	size_t dim = lhs.size();
	U s = U(0);
	for (size_t i = 0; i < dim; ++i)
		s += lhs(i) * rhs(i);
	return s;
}

template<class U>
Vector<U> operator*(double s, const Vector<U>& v) {
	size_t dim = v.size();
	Vector<U> w(dim);
	for (size_t i = 0; i < dim; ++i)
		w(i) = v(i) * s;
	return w;
}

template<class U>
Vector<U> operator*(const Vector<U>& v, double s) {
	size_t dim = v.size();
	Vector<U> w(dim);
	for (size_t i = 0; i < dim; ++i)
		w(i) = v(i) * s;
	return w;
}


/******************************************************************************/
/* The vector class declaration 																*/
/******************************************************************************/
template<class T>
class Vector {
public:

Vector() : dim(0) {
   ;
}

Vector(size_t n) : dim(n) {
	x = new T[n];
}


Vector(const Vector<T>& A) {
	dim = A.size();
	if (dim > 0) x = new T[dim];
	for (size_t i = 0; i < dim; ++i) *(x + i) = A(i);
}

Vector<T>& operator=(const Vector<T>& A) {
	dim = A.size();
	if (dim > 0) x = new T[dim];
	for (size_t i = 0; i < dim; ++i) *(x + i) = A(i);
	return *this;
}

Vector(size_t n, const T v) : dim(n) {
	x = new T[n];
	for (size_t i = 0; i < dim; ++i) *(x + i) = v;
}


size_t size() const {
	return dim;
}

void fill(T value) {
	for (size_t i = 0; i < dim; ++i) *(x + i) = value;
}

void resize(size_t n) {
	if (dim > 0) delete[] x;
	dim = n;
	x = new T[n];
}

T operator()(size_t i) const {
	return *(x + i);
}

T& operator()(size_t i) {
	return *(x + i);
}


Vector<T> operator+(const Vector<T>& y) const {
	Vector<T> v(dim);
	for (size_t i = 0; i < dim; ++i)
		v(i) = this->operator()(i) + y(i);
	return v;
}

Vector<T> operator-(const Vector<T>& y) const {
	Vector<T> v(dim);
	for (size_t i = 0; i < dim; ++i)
		v(i) = this->operator()(i) - y(i);
	return v;
}

/* multiplication operators */
friend T operator*<T>(const Vector& lhs, const Vector& rhs); // dot product
friend Vector operator*<T>(double s, const Vector& v); // scalar * vector
friend Vector operator*<T>(const Vector& v, double s); // vector * scalar

Vector<T> operator/(const T& s) const {
	Vector<T> v(dim);
	for (size_t i = 0; i < dim; ++i)
		v(i) = this->operator()(i) / s;
	return v;
}

/* returns a block of the vector*/
Vector<T> block(size_t i_end) const {
	Vector<T> v(i_end + 1);
	for (size_t i = 0; i <= i_end; ++i)
		v(i) = *(x + i);
	return v;
}

Vector<T> block(size_t i_start, size_t i_end) const {
	Vector<T> v(i_end - i_start + 1);
	for (size_t i = 0; i <= i_end - i_start; ++i)
		v(i) = *(x + i_start + i);
	return v;
}

/* return the euclidean norm of the vector */
double norm() const {
	return std::sqrt((*this)*(*this));
}

/* change the position of two elements in the vector */
void swapElements(size_t i_1, size_t i_2) {
	T e_1 = *(x + i_1);
	*(x + i_1) = *(x + i_2);
	*(x + i_2) = e_1;
}

/* returns the biggest element (absolute value) of the vector */
T biggestAbs() const {
	double big = std::abs(*x);
	for (size_t i = 1; i < dim; ++i)
		if (std::abs(*(x + i)) > big) big = std::abs(*(x + i));
	return big;
}

/* returns the index of the biggest element (absolute value) of the vector */
size_t biggestAbsIndex() const {
	size_t i = 0;
	for (size_t j = 1; j < dim; ++j)
		if (std::abs(*(x + j)) > std::abs(*(x + i))) i = j;
	return i;
}

/* returns the smallest element (absolute value) of the vector */
T smallestAbs() const {
	double small = std::abs(*x);
	for (size_t i = 1; i < dim; ++i)
		if (std::abs(*(x + i)) < small) small = std::abs(*(x + i));
	return small;
}

/* returns the index of the smallest element (absolute value) of the vector */
size_t smallestAbsIndex() const {
	size_t i = 0;
	for (size_t j = 1; j < dim; ++j)
		if (std::abs(*(x + j)) < std::abs(*(x + i))) i = j;
	return i;
}

/* returns the value of the last component */
T last_v() const {
	return *(x + dim - 1);
}

/* destructor */
virtual ~Vector() {
	//std::cout << "Exterminate Vector!" << std::endl;
	if (dim > 0) delete[] x;
}

const T* begin() const {
	return x;
}

const T* end() const {
	return x+dim;
}

private:
	T* x;
	size_t dim;
};



}

#endif

