/* 
 * File:   PotentialFlowElements.hpp
 * Author: kenobi
 *
 * Created on July 17, 2014, 11:09 AM
 */


#ifndef MATRIX_HPP
#define MATRIX_HPP

#include"Vector.hpp"

#include<iostream>
#include<cmath>

namespace nde {

    template<class T>
    class Matrix {
    public:

        /* constructors */
        Matrix() : rows(0), cols(0) {
            ;
        }

        Matrix(int numrows, int numcols) : rows(numrows), cols(numcols) {
            mat.resize(rows);
            for (int i = 0; i < rows; ++i)
                mat(i) = new Vector<T > (cols);
        }

        Matrix(const Matrix<T>& A) {
            rows = A.rows;
            cols = A.cols;
            mat.resize(rows);
            for (int i = 0; i < rows; ++i) {
                mat(i) = new Vector<T > (cols);
                for (int j = 0; j < cols; ++j)
                    mat(i)->operator()(j) = A.mat(i)->operator()(j);
            }
        }

        void resize(int numrows, int numcols) {
            rows = numrows;
            cols = numcols;
            mat.resize(rows);
            for (int i = 0; i < rows; ++i)
                mat(i) = new Vector<T > (cols);
        }

        Matrix<T>& operator=(const Matrix<T>& A) {
            rows = A.rows;
            cols = A.cols;
            mat.resize(rows);
            for (int i = 0; i < rows; ++i) {
                mat(i) = new Vector<T > (cols);
                for (int j = 0; j < cols; ++j)
                    mat(i)->operator()(j) = A(i,j);
            }
            return *this;
        }

        /* destructor */
        virtual ~Matrix() {
            //std::cout << "Exterminate Matrix!" << std::endl;
            for (int i = 0; i < rows; ++i)
                delete mat(i);
        }

        /* fill everything with a constant value */
        void fill(const T& value) {
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j)
                    this->operator()(i, j) = value;
            }
        }

        /* access an element */
        T operator()(int i, int j) const {
            return mat(i)->operator()(j);
        }

        /* modify an element */
        T& operator()(int i, int j) {
            return mat(i)->operator()(j);
        }

        /* access a row */
        Vector<T> row(int i) const {
            return *mat(i);
        }

        /* access a column */
        Vector<T> col(int i) const {
            Vector<T> x(this->numrows());
            for (int j = 0; j<this->numrows(); ++j)
                x(j) = mat(j)->operator()(i);
            return x;
        }

        /* swap two rows */
        void swapRows(int i_1, int i_2) {
            Vector<T> v_1 = *mat(i_1);
            *mat(i_1) = *mat(i_2);
            *mat(i_2) = v_1;
        }

        /* return number rows */
        int numrows() const {
            return rows;
        }

        /* return number columns */
        int numcols() const {
            return cols;
        }

        /* multiply by a scalar */
        Matrix<T> operator*(const T& s) const {
            Matrix<T> x(this->numrows(), this->numcols());
            for (int i = 0; i < this->numrows(); ++i) {
                for (int j = 0; j < this->numcols(); ++j)
                    x(i, j) = this->operator()(i, j) * s;
            }
            return x;
        }

        /* multiply by a vector */
        Vector<T> operator*(const Vector<T>& s) const {
            Vector<T> x(this->numrows());
            for (int i = 0; i < this->numrows(); ++i)
                x(i) = *mat(i) * s;
            return x;
        }

        /* multiply by a matrix */
        Matrix<T> operator*(const Matrix<T>& that) const {
            Matrix<T> x(this->numrows(), that->numcols());
            for (int i = 0; i<this->numrows(); ++i) {
                for (int j = 0; j < that->numcols(); ++j)
                    x(i, j) = this->row(i) * that.col(j);
            }
            return x;
        }

        /* solve a linear system of equations 
           A x = b */
        Vector<T> solve(const Vector<T>& b) const {

            /* compute the LU Matrix */

            /* initialize LU Matrix */
            Matrix<T> LU = *this;

            /* row scaling: largest element in row */
            Vector<T> big(rows);
            for (int i = 0; i < rows; ++i)
                big(i) = mat(i)->biggestAbs();

            /* initialize permutation vector */
            Vector<int> p(rows);
            for (int i = 0; i < rows; ++i)
                p(i) = i;

            for (int j = 0; j < rows; ++j) { /* loop through columns */

                /* upper triangular elements */
                for (int i = 1; i <= j; ++i)
                    LU(i, j) -= LU.row(i).block(i - 1) * LU.col(j).block(i - 1);

                /* lower triangular elements */
                if (j > 0) {
                    for (int i = j + 1; i < rows; ++i)
                        LU(i, j) -= LU.row(i).block(j - 1) * LU.col(j).block(j - 1);
                }

                /* row permutation by the biggest one */
                Vector<T> pcol = LU.col(j).block(j, rows - 1);
                Vector<T> pbig = big.block(j, rows - 1);
                Vector<T> merit(rows - j);
                for (int i = 0; i < rows - j; ++i)
                    merit(i) = std::abs(pcol(i)) / std::abs(pbig(i));

                int biggest_row_index = j + merit.biggestAbsIndex();

                if (biggest_row_index != j) {
                    LU.swapRows(biggest_row_index, j);
                    p.swapElements(biggest_row_index, j);
                }

                /* scale the column */
                for (int i = j + 1; i < rows; ++i)
                    LU(i, j) /= LU(j, j);

            }

            /* forward substitution */
            Vector<T> y(rows);
            y(0) = b(p(0));
            for (int i = 1; i < rows; ++i)
                y(i) = b(p(i)) - LU.row(i).block(i - 1) * y.block(i - 1);

            /* backwards substitution and solution */
            Vector<T> x(rows);
            x(rows - 1) = y(rows - 1) / LU(rows - 1, rows - 1);
            for (int i = rows - 2; i >= 0; --i)
                x(i) = (y(i) - LU.row(i).block(i + 1, rows - 1) * x.block(i + 1, rows - 1)) / LU(i, i);

            return x;

        }

    private:
        Vector<Vector<T>* > mat;
        int rows;
        int cols;
    };

}

#endif
