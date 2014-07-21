#include"combinatorial.hpp"

/*****************************************************************************/
/*		nchoosek
/*****************************************************************************/
unsigned int Combinatorial::nchoosek(unsigned int n, unsigned int k)
{

	if (k >  n) return 0;
	
	if (k*2 > n) k = n-k;

    if (k == 0) return 1;

    unsigned int result = n;
    for( unsigned int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;

}

/*****************************************************************************/
/*		nchoosek
/*****************************************************************************/
matrix<unsigned int> Combinatorial::nchoosek(vector<unsigned int> v, unsigned int k)
{


	unsigned int n = v.size();

	if (k == 0)		 {
		throw "choose value of zero not allowed-in Combinatorial::nchoosek";
	}
	else if (k > n)	 {
		throw "choose value greater than number of total elements \
				 not allowed-in Combinatorial::nchoosek";
	}
	else if (k == n) {

		matrix<unsigned int> A(1,n);
		//A.row(0) = v;
		return A;
	}
	else if (k == 1) {

		matrix<unsigned int> A(n,1);
		//A.col(0) = v;
		return A;
	}
	else  {

		unsigned int dim = nchoosek(n,k); //total number of combinations
		matrix<unsigned int> A(dim,k);  //result matrix
	
		// unsigned int cont = 0;
		// for (unsigned int l=0; l<=n-k; ++l)	{
		
			// vector<unsigned int> v1 = v.tail(n-l-1);
			// matrix<unsigned int> A1 = nchoosek(v1,k-1);

			// for (unsigned int j = 0; j<A1.norows(); ++j)	{
				// A.row(cont)[0]        = v(l);
				// A.row(cont).tail(k-1) = A1.row(j);
				// ++cont;
			// }

		// }

		return A;

	}

}

/*****************************************************************************/
/*		factorial
/*****************************************************************************/
unsigned int Combinatorial::factorial(unsigned int n)
{
	if (n==1 || n==0) return 1;
	return n*factorial(n-1);
}
