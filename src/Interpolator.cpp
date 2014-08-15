#include"curveinterpolation.hpp"
#include"combinatorial.hpp" // for Combinatorial::factorial
#include<iterator>
#include<cmath>
#include<exception> // for throw

/***************************************************************************************/
/**
 * CURVE INTERPOLATION, BASE CLASS
 * The interpolation type can be controlled by the parameter interptype.
 */
/***************************************************************************************/


/**
 * CONSTRUCTORS
 */

using namespace std;
 
CurveInterpolator::CurveInterpolator()
{ ; }

CurveInterpolator::CurveInterpolator(const vector<double>& x,
									 const vector<double>& y,
									 InterpType interptype)
									:_xvalues(x),
									 _yvalues(y),
									 _interptype(interptype)
{

	// number of nodes
	unsigned int n = _xvalues.size();
	//QL_ENSURE(y.size()==n,"Not the same number of ordinates and abscissas values");

	// calculate the spline coefficients according to the type of interpolation
	switch(_interptype)	{

	case PIECEWISECONSTANT:
	{
		_splinecoeffs = setupspline_pwc(_xvalues, _yvalues);
		_yvalues = setupordinates_pwc(_yvalues);
		break;
	}
	case LINEAR:
		_splinecoeffs = setupspline_linear(_xvalues, _yvalues);
		break;
	case CATROM:
		_splinecoeffs = setupspline_catmull(_xvalues, _yvalues);
		break;
	case SPLINEMONOTONE:
		_splinecoeffs = setupspline_monotonespline(_xvalues, _yvalues);
		break;
	default:
		throw "CurveInterpolator type not allowed";

	}

	// set up flat extrapolation for tail values by default.
	_extrapvalues.resize(2);
	_extrapvalues[0] = _yvalues[0];
	_extrapvalues[1] = _yvalues[n-1];

}

CurveInterpolator::CurveInterpolator(const vector<double>& x,
								     const matrix<double>& splinecoeffs,
									 const vector<double>& extrapvalues)
								   : _xvalues(x),
									 _splinecoeffs(splinecoeffs),
									 _extrapvalues(extrapvalues)
{

	// number of function nodes
	unsigned int n = x.size();

	// check for proper number of spline coefficients: there should be one set per interval!
	//QL_ENSURE(_splinecoeffs.nocolumns()==n-1, "bad number of spline coefficients");

	// recalculate the ordinates values
	// and set the values for the new (x,y) pairs.

	_yvalues.resize(n);
	_yvalues[0] = 0.5*(calcsplinevalue(0,_xvalues[0])+_extrapvalues[0]);
	for (unsigned int i=0; i<n-2; ++i) {
		double yl = calcsplinevalue(i,_xvalues[i+1]);
		double yr = calcsplinevalue(i+1,_xvalues[i+1]);
		_yvalues[i+1] = 0.5*(yl+yr);
	}
	_yvalues[n-1] = 0.5*(calcsplinevalue(n-2,_xvalues[n-1])+_extrapvalues[1]);

}

/**
 * Interpolates the spline value for the a given abscissas "x"
 */
double CurveInterpolator::calcvalue(double x) const
{

	// This will be the returned value.
	// Initialise it to zero.
	double y = 0.0;

	if (x<_xvalues[0])
		y = _extrapvalues[0];
	else if (x==_xvalues[0])
		y = _yvalues[0];
	else if (x>_xvalues[_xvalues.size()-1])
		y = _extrapvalues[1];
	else { // if it is within the bounds, interpolate with the polynomial.

		vector<double>::const_iterator itx = lower_bound(_xvalues.begin(), _xvalues.end(), x);
		itx--;

		const unsigned int index  = itx - _xvalues.begin();

	   //Important note!: spline values could be different from
	   //node values because of discontinuity of the spline at the nodes.
	   //When the spline is discontinuous, the function value is one half
	   //the value from the left and the value from the right.

		if (x==_xvalues[index])
			y = _yvalues[index];
		else if (x==_xvalues[index+1])
			y = _yvalues[index+1];
		else
			y = calcsplinevalue(index,x);

	}

	return y;

}



double CurveInterpolator::calcsplinevalue(unsigned int interval, double x) const
{
	double dx = x - _xvalues[interval];
	double t  = 1.0;
	double y  = 0.0;
	for (unsigned int j=0; j<_splinecoeffs.norows(); ++j) {
		y += _splinecoeffs(j,interval)*t;
		t *= dx;
	}
	return y;
}


/**
 * Interpolates the spline value for an array of abscissas values "x".
 */
vector<double> CurveInterpolator::calcvalue_v(const vector<double> &x) const
{
     vector<double> y(x.size());
     for (unsigned int i=0; i<x.size(); ++i)
           y[i] = calcvalue(x[i]);
     return y;
}

/**
 * Calculate the first derivative of the spline for a value "x" of the abscissas.
 */
double CurveInterpolator::calcdervalue(double x) const
{

	// This will be the returned value of the derivative.
	// Initialised to zero.
	double y = 0.0;


	if (x<_xvalues[0] || x>_xvalues[_xvalues.size()-1])
		y = 0.0; // The slope will be zero because of the flat extrapolation.
	else {

		unsigned int index;
		if (x==_xvalues[0])
			index = 0;
		else {
			vector<double>::const_iterator itx = lower_bound(_xvalues.begin(), _xvalues.end(), x);
			itx--;
			index  = itx - _xvalues.begin();
		}

		const unsigned int maxorder = _splinecoeffs.norows();

		// calculate the first derivative.
		double dx = x - _xvalues[index];
		double t = 1.0;
		y = 0.0;
		for (unsigned int i=1; i<maxorder; ++i) {
			y += _splinecoeffs(i,index)*(double)i*t;/*notice row major orientation
													  in _splinecoeffs*/
			t *= dx;
		}

		// Check for cases when the interpolation point lies on one of the function nodes.
		// For those cases, always return the averaged value of the first derivative approaching
		// from the left and that from the right.

		// The already calculated value above in the code is the value approaching from the left
		// if it coincides with the beginning of the interval of the value from the right if it
		// coincides with the end of the interval.

		if (x==_xvalues[index]) { //beginning of the interval -> average it with the value of
								  //the derivative at the end of the previous interval
			if (index==0)
				y = 0.5*y; //because of the flat extrapolation the value of the derivative
						   //from the right is zero.
			else {
				dx = x - _xvalues[index-1];
				t = 1.0;
				double yr = 0.0;
				for (unsigned int i=1; i<maxorder; ++i) {
					yr += _splinecoeffs(i,index-1)*(double)i*t;//notice row major orientation
																//in _splinecoeffs
					t *= dx;
				}
				y = 0.5*(y+yr);
			}
		} else if (x==_xvalues[index+1]) { //end of the interval -> average it with the value of
										   //the derivative at the beginning of the next interval
			if (index==_xvalues.size()-2)
				y = 0.5*y; //because of the flat extrapolation the value of the derivative is zero.
			else if (maxorder>1) {
				double yl = _splinecoeffs(1,index+1);
				y = 0.5*(y+yl);
			}
		}

	}

	return y;

}


/**
 * Calculate the first derivative of the spline for an array of values "x" of the abscissas.
 */
vector<double> CurveInterpolator::calcdervalue_v(const vector<double> &x) const
{
     vector<double> y(x.size());
     for (unsigned int i=0; i<x.size(); ++i)
           y[i] = calcdervalue(x[i]);
     return y;
}


/**
 * Sensitivity to the spline node values "y_i" for the node ordinates values.
 */
vector<double> CurveInterpolator::calcdervalue_yinterpnodes(double x) const
{

	//QL_ENSURE(_interptype==InterpType::LINEAR, "Method not implemented for interpolation other than linear");

	unsigned int n = _xvalues.size(); //auxiliary size variable
	vector<double> result(n, 0.0);

	if (x>=_xvalues[n-1])
		result[n-1] = 1.0;
	else if (x<=_xvalues[0])
		result[0]   = 1.0;
	else {
		vector<double >::const_iterator itx = lower_bound(_xvalues.begin(), _xvalues.end(), x);
		itx--;
		const unsigned int index  = itx - _xvalues.begin();//x is between index and index+1

		const double dx = _xvalues[index+1]-_xvalues[index];
		result[index]   = (_xvalues[index+1] - x)/dx;
		result[index+1] = ( x  - _xvalues[index])/dx;

	}

	return result;

}


/**
 * Property GETS
 */

vector<double> CurveInterpolator::getextrapvalues() const
{ return _extrapvalues; }

matrix<double> CurveInterpolator::getsplinecoeffs() const
{ return _splinecoeffs; }


/**
 * Property SETS
 */
void CurveInterpolator::setnew_x(const vector<double> &x)
{
	//QL_ENSURE(_xvalues.size()==x.size(), "New values for x cannot be loaded");
	CurveInterpolator(x,_yvalues,_interptype);
}

void CurveInterpolator::setnew_y(const vector<double> &y)
{
	//QL_ENSURE(_yvalues.size()==y.size(), "New values for y cannot be loaded");
	CurveInterpolator(_xvalues,y,_interptype);
}

void CurveInterpolator::setnew_interptype(InterpType interptype)
{
	CurveInterpolator(_xvalues,_yvalues,interptype);
}

/**
 * Set new spline coefficients for a given interval.
 */
void CurveInterpolator::setnew_splinecoeffs(unsigned int interval, const vector<double>& new_coeffs)
{

	//QL_ENSURE(new_coeffs.size()==_splinecoeffs.norows(),"Bad number of spline coefficients");

	for (unsigned int j=0; j<_splinecoeffs.norows(); ++j)
		_splinecoeffs(j, interval) = new_coeffs[j];

	// refresh now with the new spline values.
	CurveInterpolator(_xvalues,_splinecoeffs,_extrapvalues);

}

/**
 * Set a new matrix of spline coefficients.
 * @param new_coeffs: from (0 to Number of Nodes - 2) x (0 to polynomial order)
 */
void CurveInterpolator::setnew_splinecoeffs(const matrix<double>& new_coeffs)
{
	//QL_ENSURE(new_coeffs.norows() == _splinecoeffs.norows()
	//	      && new_coeffs.nocolumns() == _splinecoeffs.nocolumns(),
	//		  "Bad number of spline coefficients");

	// refresh now with the new spline values.
	CurveInterpolator(_xvalues,new_coeffs,_extrapvalues);
}


/**
 * Functions providing UTILITIES and AUXILIAR FUNCTIONALITIES
 */

/**
	Check if a given interval is positive.
*/
bool CurveInterpolator::is_positive(unsigned int interval_index) const
{

	//QL_ENSURE((interval_index >= 0) && (interval_index < _xvalues.size()),
	//		  "Interval index outside spline range");

	double a  = _xvalues[interval_index];
	double b  = _xvalues[interval_index+1];
	double dx = b-a;

	// evaluate the spline at the extreme points in the interval
	double fa = calcsplinevalue(interval_index,a);
	double fb = calcsplinevalue(interval_index,b);

	bool pos_extreme_vals = fa>=0.0 && fb>=0.0;

	if (fa<0.0 || fb<0.0) return false;

	switch (_splinecoeffs.norows()-1) { //select case by polynomial order
	case 0:
	case 1: // constant and linear cases
		return true; // nothing to check for, otherwise it would have returned false before!
		break;
	case 2: // parabola
	{
		double m = a - 0.5*_splinecoeffs(1,interval_index)/_splinecoeffs(2,interval_index);

		if (m<b && m>a) { //there is a minimum/maximum of the parabola inside this interval
			double fm = _splinecoeffs(0,interval_index)
						- 0.25 * _splinecoeffs(1,interval_index) *_splinecoeffs(1,interval_index)
						       / _splinecoeffs(2,interval_index);
			return fm>=0.0;
		} else
			return true;

		break;
	}
	case 3: // cubic function
	{
		double r = _splinecoeffs(2,interval_index)/_splinecoeffs(1,interval_index);
		double w = r*r - 3.0*_splinecoeffs(3,interval_index)
		                    /_splinecoeffs(1,interval_index);

		if (w>=0.0) { //the cubic function has a local maximum and minimum

			double m1 = a - r + sqrt(w);
			double m2 = a - r - sqrt(w);
			double fm1 = fa, fm2 = fb; // Initialised to any positive value

			// check if maximum and/or minimum is inside this interval
			if (m1<b && m1>a) fm1 = calcvalue(m1);
			if (m2<b && m2>a) fm2 = calcvalue(m2);

			return fm1>=0.0 && fm2>=0.0;

		} else
			return true;

		break;
	}
	default:
		throw "Positive check not implemented for polynomial order greater than three";
		break;
	}

}

/**
	Check if the spline is positive all over the place.
*/
bool CurveInterpolator::is_positivespline() const
{
	bool pos_allintervals = true;
	for (unsigned int k=0; k<_xvalues.size()-1; ++k)
		pos_allintervals *= is_positive(k);
	return pos_allintervals;
}


/**
	Add two splines classes
	i.e.: this routine adds the spline coefficients
*/
CurveInterpolator CurveInterpolator::operator+(const  CurveInterpolator& other_spline) const
{

	matrix<double> other_splinecoeffs = other_spline.getsplinecoeffs();

	//QL_ENSURE(other_splinecoeffs.nocolumns()==_splinecoeffs.nocolumns(),
	//		  "the two splines cannot be added up under the current implementation \
	//		   because they do not share the same intervals.");

	unsigned int order = max(_splinecoeffs.norows(),other_splinecoeffs.norows());

	matrix<double> new_splinecoeffs(order, _splinecoeffs.nocolumns());

	for (unsigned int j=0; j<_splinecoeffs.nocolumns(); ++j)  {
		for (unsigned int k=0; k<_splinecoeffs.norows(); ++k)
			new_splinecoeffs(k,j) = _splinecoeffs(k,j);
		for (unsigned int k=0; k<other_splinecoeffs.norows(); ++k)
			new_splinecoeffs(k,j) += other_splinecoeffs(k,j);
	}


	vector<double> other_extrapvalues = other_spline.getextrapvalues();
	vector<double> new_extrapvalues(2);
	for (unsigned int j=0; j<2; ++j)
		new_extrapvalues[j] = _extrapvalues[j] + other_extrapvalues[j];

	return  CurveInterpolator(_xvalues, new_splinecoeffs, new_extrapvalues);

}

/**
	Subtract two splines classes
	i.e.: this routine subtracts the spline coefficients
*/
CurveInterpolator CurveInterpolator::operator-(const  CurveInterpolator& other_spline) const
{

	matrix<double> other_splinecoeffs = other_spline.getsplinecoeffs();

	//QL_ENSURE(other_splinecoeffs.nocolumns()==_splinecoeffs.nocolumns(),
	//		  "the two splines cannot be subtracted under the current implementation \
	//		   because they do not share the same intervals.");

	unsigned int order = max(_splinecoeffs.norows(),other_splinecoeffs.norows());

	matrix<double> new_splinecoeffs(order, _splinecoeffs.nocolumns());

	for (unsigned int j=0; j<_splinecoeffs.nocolumns(); ++j)  {
		for (unsigned int k=0; k<_splinecoeffs.norows(); ++k)
			new_splinecoeffs(k,j) = _splinecoeffs(k,j);
		for (unsigned int k=0; k<other_splinecoeffs.norows(); ++k)
			new_splinecoeffs(k,j) -= other_splinecoeffs(k,j);
	}

	vector<double> other_extrapvalues = other_spline.getextrapvalues();
	vector<double> new_extrapvalues(2);
	for (unsigned int j=0; j<2; ++j)
		new_extrapvalues[j] = _extrapvalues[j] - other_extrapvalues[j];

	return  CurveInterpolator(_xvalues, new_splinecoeffs, new_extrapvalues);

}

/**
 * Scale up a spline
 * @param factor: a double factor to multiply by
*/
CurveInterpolator CurveInterpolator::operator*(double factor) const
{

	matrix<double> new_splinecoeffs(_splinecoeffs.norows(),_splinecoeffs.nocolumns());

	for (unsigned int j=0; j<_splinecoeffs.norows(); ++j)  {
		for (unsigned int k=0; k<_splinecoeffs.nocolumns(); ++k)
			new_splinecoeffs(j,k) = factor*_splinecoeffs(j,k);
	}

	vector<double> new_extrapvalues(2);
	for (unsigned int j=0; j<2; ++j)
		new_extrapvalues[j] = _extrapvalues[j]*factor;

	return  CurveInterpolator(_xvalues, new_splinecoeffs, new_extrapvalues);

}

/**
 * Scale down a spline
 * @param factor: a double factor to divide by
 */
CurveInterpolator CurveInterpolator::operator/(double factor) const
{

	matrix<double> new_splinecoeffs(_splinecoeffs.norows(),_splinecoeffs.nocolumns());

	for (unsigned int j=0; j<_splinecoeffs.norows(); ++j)  {
		for (unsigned int k=0; k<_splinecoeffs.nocolumns(); ++k)
			new_splinecoeffs(j,k) = _splinecoeffs(j,k)/factor;
	}

	vector<double> new_extrapvalues(2);
	for (unsigned int j=0; j<2; ++j)
		new_extrapvalues[j] = _extrapvalues[j]/factor;

	return  CurveInterpolator(_xvalues, new_splinecoeffs, new_extrapvalues);

}

/**
 * Multiply two splines sharing the same intervals
 */
 CurveInterpolator CurveInterpolator::operator*(const CurveInterpolator& other_spline) const
{

	unsigned int num_intervals = _splinecoeffs.nocolumns();

	matrix<double> other_splinecoeffs = other_spline.getsplinecoeffs();
	//QL_ENSURE(other_splinecoeffs.nocolumns()==num_intervals, \
	//		  "The two splines cannot be multiplied because they don't share the same number of intervals");

	unsigned int order       = _splinecoeffs.norows()-1;
	unsigned int other_order = other_splinecoeffs.norows()-1;

	matrix<double> new_splinecoeffs(order+other_order+1,_splinecoeffs.nocolumns());
	new_splinecoeffs.fill(0.0);

	for (unsigned int l=0; l<num_intervals; ++l) {
		for (unsigned int j=0; j<=order; ++j) {
			for (unsigned int k=0; k<=other_order; ++k)
				new_splinecoeffs(j+k,l) += _splinecoeffs(j,l)*other_splinecoeffs(k,l);
		}
	}

	vector<double> other_extrapvalues = other_spline.getextrapvalues();
	vector<double> new_extrapvalues(2);
	for (unsigned int j=0; j<2; ++j)
		new_extrapvalues[j] = _extrapvalues[j]*other_extrapvalues[j];

	return  CurveInterpolator(_xvalues, new_splinecoeffs, new_extrapvalues);

}


/**
 * Take the n-th derivative order
 * @param order: the order of derivation, i.e. =0 no derivative, =1 first derivative ,
 * 											   =2 second derivative, ...
 */
CurveInterpolator CurveInterpolator::derive(unsigned int order) const
{

	int old_order = _splinecoeffs.norows()-1;
	int new_order = old_order - order; // it can be positive or negative
												 // negative order would mean that the function
												 // vanishes off and equals zero.

	matrix<double> new_splinecoeffs(max(new_order,(int)0.0)+(int)1.0, _splinecoeffs.nocolumns());
	new_splinecoeffs.fill(0.0);

	if (new_order>=0) {

		vector<double> factor(new_order+1);
		for (unsigned int k=0; k<=new_order; ++k)
			factor[k] = (double)Combinatorial::factorial(k+order)/ (double)Combinatorial::factorial(k);

		for (unsigned int j=0; j<_splinecoeffs.nocolumns(); ++j) {
			for (unsigned int k=0; k<=new_order; ++k)
				new_splinecoeffs(k,j) = factor[k]*_splinecoeffs(k+order,j);
		}

	}

	// the derivative in the outside of the rage is zero, because of the flat extrapolation
	vector<double> extrapvalues(2, 0.0);

	return  CurveInterpolator(_xvalues, new_splinecoeffs, extrapvalues);

}


/* THIS FUNCTIONALITY HAS TO BE RE-THOUGHT */

///**
// * Integrate the spline once, between starting point in the interval
// * up to x.
// * The resulting spline is defined as g(x) = int_{x_start}^x f(s) ds
// */
// CurveInterpolator  CurveInterpolator::integrate() const
//{
//
//	matrix<double> new_splinecoeffs(_splinecoeffs.norows()+1,_splinecoeffs.nocolumns());
//	new_splinecoeffs.fill(0.0);
//
//	double area = 0.0;
//	double v;
//
//	for (unsigned int j=0; j<_splinecoeffs.nocolumns(); ++j) {
//
//		double dx = _xvalues[j+1] - _xvalues[j];
//		v = 0.0;
//		double t = dx;
//		for (unsigned int k=1; k<_splinecoeffs.norows()+1; ++k) {
//			new_splinecoeffs(k,j) = _splinecoeffs(k-1,j)/k;
//			v  += new_splinecoeffs(k,j)*t;
//			t  *= dx;
//		}
//
//		new_splinecoeffs(j,0) = area;
//		area += v;
//
//	}
//
//	return  CurveInterpolator(_xvalues,new_splinecoeffs);
//
//}





/***************************************************************************************/
/*
	Auxiliary 1D interpolation functions
*/
/***************************************************************************************/




/*
	Calculate interpolation coefficients for Piece-Wise Constant
*/
matrix<double> setupspline_pwc(const vector<double>& x, const vector<double>& y)
{

	const unsigned int n = x.size();
	matrix<double> interpcoeffs(1,n-1);

	for (unsigned int i=0; i<n-1; ++i)
		interpcoeffs(0,i) = y[i];

	return interpcoeffs;

}


/**
 * Adjust the PWC node ordinates so that they match the criteria that nodes values
 * are the averaged values of the spline from the left and the right sides of the interval
 */
vector<double> setupordinates_pwc(const vector<double>& y)
{
	vector<double> y_avg(y.size());

	y_avg[0] = y.front();
	for (unsigned int j=1; j<y.size(); ++j)
		y_avg[j] = 0.5*(y[j-1]+y[j]);
	y_avg[y.size()-1] = y.back();

	return y_avg;
}

/*
	Calculate interpolation coefficients for Linear
*/
matrix<double> setupspline_linear(const vector<double>& x, const vector<double>& y)
{

	const unsigned int n = x.size();
	matrix<double> interpcoeffs(2,n-1);

	for (unsigned int i=0; i<n-1; ++i)	{
		interpcoeffs(0,i) = y[i];
		interpcoeffs(1,i) = (y[i+1]-y[i])/(x[i+1] -x[i]);
	}

	return interpcoeffs;

}


/*
	Calculate interpolation coefficients for CATMUL-ROM
*/
matrix<double> setupspline_catmull(const vector<double>& x, const vector<double>& y)
{

	const unsigned int n = x.size();
	matrix<double> interpcoeffs(4,n-1);

	for (unsigned int i=0; i<n-1; ++i) {

		double l  = x[i+1] - x[i];
		double l2 = l*l;
		double l3 = l*l2;
		double sa;
		double sb;

		if (i>0)
			sa  = (y[i+1] -y[i-1]);
		else
			sa  = (y[1]-y[0]);

		if (i<n-2)
			sb  = (y[i+2]-y[i]);
		else
			sb  = (y[n-1]-y[n-2]);

		// now go with the interpolation coefficients.

		// order zero
		interpcoeffs(0,i) = y[i];

		// order one
		interpcoeffs(1,i) = sa/l;

		// order square
		interpcoeffs(2,i) = (3.0*(y[i+1] - y[i]) - 2.0*sa - sb)/l2;

		// order cubic
		interpcoeffs(3,i) = (-2.0*(y[i+1] - y[i]) + sa + sb)/l3;

	}

	return interpcoeffs;


}

/*
	Calculate interpolation coefficients for Monotone splines
*/
matrix<double> setupspline_monotonespline(const vector<double>& x, const vector<double>& y)
{

	// these are the interpolation coefficients
	// y(x, x_i-1<x<x_i) = _ai (x - x_i)^3 + _bi (x - x_i)^2 + _ci (x - x_i) + _di
	const unsigned int n = x.size();
	matrix<double> interpcoeffs(4,n-1);

	// initialize x-grid step sizes and slopes.
	vector<double> h(n-1);
	vector<double> s(n-1);
	for (unsigned int i = 0; i < n-1; ++i) {
		h[i] = x[i+1] - x[i];
		s[i] = (y[i+1] - y[i])/h[i];
	}

	if (n>2) {

		vector<double> p(n);
		p[0] = 0.0; p[n-1] = 0.0;
		for (unsigned int i = 1; i < n-1; ++i)
			p[i] = (s[i-1]*h[i] + s[i]*h[i-1])/(h[i-1] + h[i]);

		vector<double> yp(n);

		for (unsigned int i = 1; i < n-1; ++i) {

			double s1 = abs(s[i-1]);
			double s2 = abs(s[i]  );

			if (s[i-1]*s[i] <= 0.0) {
				yp[i] = 0.0;
			}else if ( ( abs(p[i]) > 2.0*s1 ) ||
						( abs(p[i]) > 2.0*s2 ) ) {
				double a = (s1 > 0.0) - (s1 < 0.0); //sign function.
				yp[i] = 2.0*a*min(s1,s2);
			}else {
				yp[i] = p[i];
			}

		}

		// second derivative at extreme points equal to zero boundary condition.
		yp[0]   = 1.5*s[0]   - 0.5*yp[1];
		yp[n-1] = 1.5*s[n-2] - 0.5*yp[n-3];

		for (unsigned int i = 0; i < n-1; ++i) {
			interpcoeffs(3,i) = (yp[i] + yp[i+1] - 2.0*s[i])/(h[i]*h[i]);
			interpcoeffs(2,i) = (3.0*s[i] - 2.0*yp[i] - yp[i+1])/h[i];
			interpcoeffs(1,i) = yp[i];
			interpcoeffs(0,i) = y[i];
		}

	}
	else {

		// in this case, with second derivative at extreme points equal to zero,
		// the solution will be just the linear interpolation.

		interpcoeffs(3,0) = 0.0; interpcoeffs(3,1) = 0.0;
		interpcoeffs(2,0) = 0.0; interpcoeffs(2,1) = 0.0;
		interpcoeffs(1,0) = s[0]; interpcoeffs(1,1) = s[1];
		interpcoeffs(0,0) = y[0]; interpcoeffs(0,1) = y[1];

	}

	return interpcoeffs;

}
