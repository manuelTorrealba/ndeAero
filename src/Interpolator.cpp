#include"Interpolator.hpp"

#include<algorithm>
//#include<iterator>
#include<cmath>
#include<exception> // for throw


using namespace std;

namespace nde {

    Interpolator1D::Interpolator1D(const Vector<double>& x,
            const Vector<double>& y,
            InterpolationType interp_type)
    : x_values(x),
    y_values(y),
    interpolation_type(interp_type) {

        // number of nodes
        unsigned int n = x_values.size();

        // calculate the spline coefficients according to the type of interpolation
        switch (interpolation_type) {

            case LINEAR:
                spline_coeffs = setupSplineLinear(x_values, y_values);
                break;
            case SPLINE_MONOTONE:
                spline_coeffs = setupSplineMonotoneSpline(x_values, y_values);
                break;
            default:
                throw "CurveInterpolator type not allowed";

        }

        // set up flat extrapolation for tail values by default.
        extrap_values.resize(2);
        extrap_values(0) = y_values(0);
        extrap_values(1) = y_values(n - 1);

    }

    /**
     * Interpolates the spline value for the a given abscissas "x"
     */
    double Interpolator1D::calcValue(double x) const {

        // This will be the returned value.
        // Initialize it to zero.
        double y = 0.0;

        if (x < x_values(0))
            y = extrap_values(0);
        else if (x == x_values(0))
            y = y_values(0);
        else if (x > x_values(x_values.size() - 1))
            y = extrap_values(1);
        else { // if it is within the bounds, interpolate with the polynomial.

            const double* itx = lower_bound(x_values.begin(), x_values.end(), x);
            itx--;

            const unsigned int index = itx - x_values.begin();

            //Important note!: spline values could be different from
            //node values because of discontinuity of the spline at the nodes.
            //When the spline is discontinuous, the function value is one half
            //the value from the left and the value from the right.

            if (x == x_values(index))
                y = y_values(index);
            else if (x == x_values(index + 1))
                y = y_values(index + 1);
            else
                y = calcSplineValue(index, x);

        }

        return y;

    }

    double Interpolator1D::calcSplineValue(int interval, double x) const {
        double dx = x - x_values(interval);
        double t = 1.0;
        double y = 0.0;
        for (int j = 0; j < spline_coeffs.numrows(); ++j) {
            y += spline_coeffs(j, interval) * t;
            t *= dx;
        }
        return y;
    }




    /***************************************************************************************/
    /*
            Auxiliary 1D interpolation functions
     */
    /***************************************************************************************/

    /*
            Calculate interpolation coefficients for Linear
     */
    Matrix<double> setupSplineLinear(const Vector<double>& x, const Vector<double>& y) {

        const unsigned int n = x.size();
        Matrix<double> interp_coeffs(2, n - 1);

        for (unsigned int i = 0; i < n - 1; ++i) {
            interp_coeffs(0, i) = y(i);
            interp_coeffs(1, i) = (y(i + 1) - y(i)) / (x(i + 1) - x(i));
        }

        return interp_coeffs;

    }

    /*
            Calculate interpolation coefficients for Monotone splines
     */
    Matrix<double> setupSplineMonotoneSpline(const Vector<double>& x, const Vector<double>& y) {

        // these are the interpolation coefficients
        // y(x, x_i-1<x<x_i) = _ai (x - x_i)^3 + _bi (x - x_i)^2 + _ci (x - x_i) + _di
        const unsigned int n = x.size();
        Matrix<double> interp_coeffs(4, n - 1);

        // initialize x-grid step sizes and slopes.
        Vector<double> h(n - 1);
        Vector<double> s(n - 1);
        for (unsigned int i = 0; i < n - 1; ++i) {
            h(i) = x(i + 1) - x(i);
            s(i) = (y(i + 1) - y(i)) / h(i);
        }

        if (n > 2) {

            Vector<double> p(n);
            p(0) = 0.0;
            p(n - 1) = 0.0;
            for (unsigned int i = 1; i < n - 1; ++i)
                p(i) = (s(i - 1) * h(i) + s(i) * h(i - 1)) / (h(i - 1) + h(i));

            Vector<double> yp(n);

            for (unsigned int i = 1; i < n - 1; ++i) {

                double s1 = abs(s(i - 1));
                double s2 = abs(s(i));

                if (s(i - 1) * s(i) <= 0.0) {
                    yp(i) = 0.0;
                } else if ((abs(p(i)) > 2.0 * s1) ||
                        (abs(p(i)) > 2.0 * s2)) {
                    double a = (s1 > 0.0) - (s1 < 0.0); //sign function.
                    yp(i) = 2.0 * a * min(s1, s2);
                } else {
                    yp(i) = p(i);
                }

            }

            // second derivative at extreme points equal to zero boundary condition.
            yp(0) = 1.5 * s(0) - 0.5 * yp(1);
            yp(n - 1) = 1.5 * s(n - 2) - 0.5 * yp(n - 3);

            for (unsigned int i = 0; i < n - 1; ++i) {
                interp_coeffs(3, i) = (yp(i) + yp(i + 1) - 2.0 * s(i)) / (h(i) * h(i));
                interp_coeffs(2, i) = (3.0 * s(i) - 2.0 * yp(i) - yp(i + 1)) / h(i);
                interp_coeffs(1, i) = yp(i);
                interp_coeffs(0, i) = y(i);
            }

        } else {

            // in this case, with second derivative at extreme points equal to zero,
            // the solution will be just the linear interpolation.

            interp_coeffs(3, 0) = 0.0;
            interp_coeffs(3, 1) = 0.0;
            interp_coeffs(2, 0) = 0.0;
            interp_coeffs(2, 1) = 0.0;
            interp_coeffs(1, 0) = s(0);
            interp_coeffs(1, 1) = s(1);
            interp_coeffs(0, 0) = y(0);
            interp_coeffs(0, 1) = y(1);

        }

        return interp_coeffs;

    }

}
