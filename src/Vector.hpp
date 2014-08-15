/* 
 * File:   PotentialFlowElements.hpp
 * Author: kenobi
 *
 * Created on July 17, 2014, 11:09 AM
 */

#ifndef VECTOR_HPP
#define VECTOR_HPP

#include<cstdlib>
#include<iostream>

namespace nde {

    template<class T>
    class Vector {
    public:

        Vector() : dim(0) {
            ;
        }

        Vector(int n) : dim(n) {
            x = (T*) std::malloc(sizeof (T) * n);
        }

        Vector(const Vector<T>& A) {
            dim = A.size();
            if (dim > 0) x = (T*) std::malloc(sizeof (T) * dim);
            for (int i = 0; i < dim; ++i) *(x + i) = A(i);
        }

        Vector<T>& operator=(const Vector<T>& A) {
            dim = A.size();
            if (dim > 0) x = (T*) std::malloc(sizeof (T) * dim);
            for (int i = 0; i < dim; ++i) *(x + i) = A(i);
            return *this;
        }

        int size() const {
            return dim;
        }

        void fill(T value) {
            for (int i = 0; i < dim; ++i) *(x + i) = value;
        }

        void resize(int n) {
            if (dim > 0) std::free(x);
            dim = n;
            x = (T*) std::malloc(sizeof (T) * n);
        }

        T operator()(int i) const {
            return *(x + i);
        }

        T& operator()(int i) {
            return *(x + i);
        }

        T operator*(const Vector<T>& y) const {
            T s = T(0);
            for (int i = 0; i < dim; ++i)
                s += *(x + i) * y(i);
            return s;
        }

        Vector<T> operator+(const Vector<T>& y) const {
            Vector<T> v(dim);
            for (int i = 0; i < dim; ++i)
                v(i) = this->operator()(i) + y(i);
            return v;
        }

        Vector<T> operator-(const Vector<T>& y) const {
            Vector<T> v(dim);
            for (int i = 0; i < dim; ++i)
                v(i) = this->operator()(i) - y(i);
            return v;
        }

        Vector<T> operator*(const T& s) const {
            Vector<T> v(dim);
            for (int i = 0; i < dim; ++i)
                v(i) = this->operator()(i) * s;
            return v;
        }

        Vector<T> operator/(const T& s) const {
            Vector<T> v(dim);
            for (int i = 0; i < dim; ++i)
                v(i) = this->operator()(i) / s;
            return v;
        }

        /* returns a block of the vector*/
        Vector<T> block(int i_end) const {
            Vector<T> v(i_end + 1);
            for (int i = 0; i <= i_end; ++i)
                v(i) = *(x + i);
            return v;
        }

        Vector<T> block(int i_start, int i_end) const {
            Vector<T> v(i_end - i_start + 1);
            for (int i = 0; i <= i_end - i_start; ++i)
                v(i) = *(x + i_start + i);
            return v;
        }

        /* return the euclidean norm of the vector */
        double norm() const {
            return std::sqrt((*this)*(*this));
        }

        /* change the position of two elements in the vector */
        void swapElements(int i_1, int i_2) {
            T e_1 = *(x + i_1);
            *(x + i_1) = *(x + i_2);
            *(x + i_2) = e_1;
        }

        /* returns the biggest element (absolute value) of the vector */
        T biggestAbs() const {
            double big = std::abs(*x);
            for (int i = 1; i < dim; ++i)
                if (std::abs(*(x + i)) > big) big = std::abs(*(x + i));
            return big;
        }

        /* returns the index of the biggest element (absolute value) of the vector */
        int biggestAbsIndex() const {
            int i = 0;
            for (int j = 1; j < dim; ++j)
                if (std::abs(*(x + j)) > std::abs(*(x + i))) i = j;
            return i;
        }

        virtual ~Vector() {
            //std::cout << "Exterminate Vector!" << std::endl;
            if (dim > 0) std::free(x);
        }

        const T* begin() const {
            return x;
        }

        const T* end() const {
            return x+dim;
        }
        
    private:
        T* x;
        int dim;
    };

}

#endif
