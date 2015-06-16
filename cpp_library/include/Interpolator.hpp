#ifndef INTERPOLATOR_HPP
#define INTERPOLATOR_HPP

/**
        Interpolation class

        HARD-CODED CONVENTIONS:
        - flat extrapolation beyond boundary points.

 */

#include"Matrix.hpp"
#include"Vector.hpp"

namespace nde {

    enum InterpolationType {
        LINEAR = 1,
        SPLINE_MONOTONE = 4
    };

    class Interpolator1D {
    public:
        /**
                Construct by entering the (x,y)-node values
                @param x			: abscissas.
                @param y			: ordinates.
                @param interptype	: interpolation type.
         */
        Interpolator1D(const Vector<double>& x,
                const Vector<double>& y,
                InterpolationType interp_type);

        /**
                Interpolate main function.
                @param x	: coordinate to interpolate to.
         */
        double operator()(double x) const;

		  /**
					Calculation of a derivate with order h.
					@param h: order of the derivative.
					@param x: coordinate to interpolate to.
			*/
			double operator()(unsigned int h, double x) const;

			/**
					Get x_min and x_max interpolation range
			 */
			Vector<double> getInterpolationRange() const;

    protected:
        InterpolationType interpolation_type; /** < interpolation type	*/
        Vector<double> x_values; /** < Vector containing x abscissas values*/
        Vector<double> y_values; /** < Vector containing y ordinates values*/

    private:
        Vector<double> extrap_values; //! _extrapvalues[0] is the low tail 
			//! value and _extrapvalues[1] is the upper tail value.
        Matrix<double> spline_coeffs; //! dimensions (from 0 to polynomial 
												  //! order, from 0 to num_nodes-2)
        /**
         * Function to calculate the spline value for a given interval at
			* a given abscissas.
         * This function needs _xvalues and _splinecoeffs properly initialised.
         */
        double calcSplineValue(unsigned int interval, unsigned int h, 
										double x) const;
    };

    /***************************************************************************************/
    /*
            Auxiliary 1D interpolation functions
     */
    /***************************************************************************************/

    /*
            Calculate interpolation coefficients for Linear
     */
    Matrix<double> setupSplineLinear(const Vector<double>& x, const Vector<double>& y);

    /*
            Calculate interpolation coefficients for Monotone splines
     */
    Matrix<double> setupSplineMonotoneSpline(const Vector<double>& x, const Vector<double>& y);


}


#endif


