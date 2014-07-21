#ifndef _CURVE_INTERPOLATION_H_
#define _CURVE_INTERPOLATION_H_

/**
	Curve interpolation generic class

	It is intended that all 1D interpolation methods are implemented inside
	this class CurveInterpolator

	Interpolation type can be controlled with the class variable _interptype

	In general, function values are obtained by evaluating the polynomial

	f(x) = sum_i c_i (x-x_i)^i, from i=0, p, where p is the grade,
												   x_i is the starting node of the i-th interval.

	The coefficients c_i can be generated for linear interpolation and different kind of splines.

	HARD-CODED CONVENTIONS:

	- function values at the nodes is equal to the averaged values of the spline
	  from the left and from the right.

	- flat extrapolation below x_0 and above x_N.



*/

#include"matrix.hpp"
#include<vector>

using namespace std;

class CurveInterpolator
{
public:

	enum InterpType {
		PIECEWISECONSTANT = 1,
		LINEAR = 2,
		CATROM = 3, 
		SPLINEMONOTONE = 4
	};
	
	/**************************************************************
	 *Constructors
	 **************************************************************/
	/**
		By default empty constructor
	*/
	CurveInterpolator();

	/**
		Construct by entering the (x,y)-node values
		@param x			: abscissas.
		@param y			: ordinates.
		@param interptype	: interpolation type.
	*/
	CurveInterpolator(const vector<double>& x,
                      const vector<double>& y,
					  InterpType interptype);

	/**
		Construct by entering the x-node values, the spline coefficients and
		the extrapolating values.
		y-node values will be recalculated as averaged values from the left and from the right.
		Notice that with this configuration the spline is not necessarily continuous.

		@param x			: abscissas.
		@param y			: ordinates.
		@param interptype	: interpolation type.
	*/
	CurveInterpolator(const vector<double>& x,
			  const matrix<double>& splinecoeffs,
			  const vector<double>& extrapvalues);

	/**************************************************************
	 *Main Functionalities
	 **************************************************************/

	/**
		Interpolate main function.
		@param x	: coordinate to interpolate to.
	*/
	double calcvalue(double x) const;

	/**
	 * Interpolate in a vector function.
	 * Implemented through looping the previous function.
	 * @param x : vector of coordinates to interpolate to.
	*/
	vector<double> calcvalue_v(const vector<double>& x) const;

	/**
	 * Calculate first derivative of the main interpolating function.
	 * @param x	: coordinate to interpolate to.
	*/
	double calcdervalue(double x) const;

	/**
	 * Calculate first derivative of the main interpolating function.
	 * Interpolate in a vector function.
	 * Implemented through looping the previous function.
	 * @param x	: vector of coordinates to interpolate to.
	*/
	vector<double> calcdervalue_v(const vector<double>& x) const;

	/**
		Calc derivative at point x with respect to node y values.
		Currently only LINEAR interpolation is supported.
		@param x	: coordinate for derivative evaluation.
	*/
	vector<double> calcdervalue_yinterpnodes(double x) const;



	/**************************************************************
	 *Property gets
	 **************************************************************/

	/**
	 * Get the extrapolating flat values.
	 */
	vector<double> getextrapvalues() const;

	/**
	 * Get the spline coefficients.
	 */
	matrix<double> getsplinecoeffs() const;


	/**************************************************************
	 *Property sets
	 **************************************************************/

	/**
		Change the coordinate x values.
		@param x	: vector of new abscissas.
	*/
	void setnew_x(const vector<double>& x);

	/**
		Change the coordinate y values.
		@param y	: vector of new ordinates.
	*/
	void setnew_y(const vector<double>& y);

	/**
		Change the interpolation type.
		@param interptype	: new interpolation method.
	*/
	void setnew_interptype(InterpType interptype);

	/**
	 * Set spline coefficients for a given interval.
	 * @param interval: from 0 to Number of Nodes - 2
	 * @param new_coeffs: from 0 to polynomial order.
	 */
	void setnew_splinecoeffs(unsigned int interval, const vector<double>& new_coeffs);

	/**
	 * Set a matrix of spline coefficients.
	 * @param new_coeffs: from (0 to Number of Nodes - 2) x (0 to polynomial order)
	 */
	void setnew_splinecoeffs(const matrix<double>& splinecoeffs);

	/**************************************************************
	 * Additional functionality
	 **************************************************************/

	/**
	 * Check if a given interval is always a positive function.
	 * @param interval_index : ranges from 0 to Number of Points - 2.
	 */
	bool is_positive(unsigned int interval_index) const;

	/**
		Check if the spline is always a positive function.
	*/
	bool is_positivespline() const;

	/**
		Add two splines classes
		i.e.: this routine adds the spline coefficients
	*/
	CurveInterpolator operator+(const CurveInterpolator& other_spline) const;

	/**
		Subtract two splines classes
		i.e.: this routine subtracts the spline coefficients
	*/
	CurveInterpolator operator-(const CurveInterpolator& other_spline) const;

	/**
	 * Scale up a spline
	 * @param factor: a double factor to multiply by
	*/
	CurveInterpolator operator*(double factor) const;

	/**
	 * Scale down a spline
	 * @param factor: a double factor to divide by
	 */
	CurveInterpolator operator/(double factor) const;

	/**
	 * Multiply two splines sharing the same intervals
	 */
	CurveInterpolator operator*(const CurveInterpolator& other_spline) const;

	/**
	 * Take the n-th derivative order
	 * @param order: the order of derivation, i.e. =0 no derivative, =1 first derivative ,
	 * 											   =2 second derivative, ...
	 */
	CurveInterpolator derive(unsigned int order) const;


	/* THIS FUNCTIONALITY HAS TO BE RETHOUGHT */
	/**
	 * Integrate the spline once, between starting point in the interval
	 * up to x.
	 * The resulting spline is defined as g(x) = int_{x_start}^x f(s) ds
	 */
	//CurveInterpolator integrate() const;


protected:
	InterpType _interptype;/** < interpolation type	*/
	vector<double> _xvalues; /** < Vector containing x abscissas values*/
	vector<double> _yvalues; /** < Vector containing y ordinates values*/

private:
	vector<double> _extrapvalues; //! _extrapvalues[0] is the low tail value and
							 //! _extrapvalues[1] is the upper tail value.
	matrix<double> _splinecoeffs; //! dimensions (from 0 to polynomial order, from 0 to num_nodes-2)

	/**
	 * Function to calculate the spline value for a given interval at a given abscissas.
	 * This function needs _xvalues and _splinecoeffs properly initialised.
	 */
	double calcsplinevalue(unsigned int interval, double x) const;

};


/**************************************************************************************/
/**
 * Some auxiliary functions to calculate spline coefficients with different
 * methodologies.
 *
 * These functions should return matrix<double> SPLINECOEFFS: (from 0 to polynomial order,
 * 														  from 0 to Number Nodes - 2)
 */
/**************************************************************************************/

/**
 * CURVE INTERPOLATION PIECE WISE CONSTANT
 */
matrix<double> setupspline_pwc(const vector<double>& xvalues, const vector<double>& yvalues);
vector<double> setupordinates_pwc(const vector<double>& yvalues);

/**
 * CURVE INTERPOLATION LINEAR
 */
matrix<double> setupspline_linear(const vector<double>& xvalues, const vector<double>& yvalues);

/**
 * CURVE INTERPOLATION WITH CATROM_MUL SPLINES
 */
matrix<double> setupspline_catmull(const vector<double>& xvalues, const vector<double>& yvalues);

/**
 * CURVE INTERPOLATION WITH MONOTONE SPLINES
 */
matrix<double> setupspline_monotonespline(const vector<double>& xvalues, const vector<double>& yvalues);


#endif
