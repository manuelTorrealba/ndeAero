/**
 * File:   TestODESolver.hpp
 * Author: kenobi
 *
 * Created on May 7, 2015, 8:53 AM
 */

#ifndef INCLUDE_TESTODESOLVER_HPP
#define INCLUDE_TESTODESOLVER_HPP

#include "ODESolver.hpp"

namespace nde {

class SimpleExponential : public ODESolver {
public:
SimpleExponential(double lambda, unsigned int order);
virtual Vector<double> odeSolverDy(double t, const Vector<double>& y) const;

private:
double _lambda;

};

class Blasius : public ODESolver {
public:
Blasius(unsigned int order);
virtual Vector<double> odeSolverDy(double t, const Vector<double>& y) const;

};


void testODESolver();


}

#endif

