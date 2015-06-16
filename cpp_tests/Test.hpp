/**
 * File:   Test.hpp
 * Author: kenobi
 *
 * Created on Jun 8, 2015, 11:00 PM
 */

#ifndef INCLUDE_TEST_HPP
#define INCLUDE_TEST_HPP

#include "Vector.hpp"

namespace nde {

	bool testThinAirfoilTheory();
	bool testODESolver();
	bool testNACAAirfoil(const Vector<unsigned int> &naca_codenum,
							const std::string& out_file_name);

}

#endif

