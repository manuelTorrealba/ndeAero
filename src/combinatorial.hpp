#include"matrix.hpp"
#include<vector>


class Combinatorial
{
public:
/**
	Number of combinations: take k elements out of N
*/
static unsigned int nchoosek(unsigned int n, unsigned int k);

/**
	Number of combinations: take k elements out of N
*/
static matrix<unsigned int> nchoosek(vector<unsigned int> v, unsigned int k);

/**
	Factorial of an integer
*/
static unsigned int factorial(unsigned int n);

};
