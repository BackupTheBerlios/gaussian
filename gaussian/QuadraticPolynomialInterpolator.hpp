/***************************************************************************
 *            QuadraticPolynomialInterpolator.hpp
 *
 *  Sun Oct  7 23:12:56 2007
 *  Copyright  2007  mjhmeyer@yahoo.com
 *  
 ****************************************************************************/

#ifndef ____QuadraticPolynomialInterpolator_H_17373___
#define ____QuadraticPolynomialInterpolator_H_17373___


#include "Typedefs.hpp"


class QuadraticPolynomialInterpolator {
	
public:
	
/** Least squares optimal interpolation \f$ Q(x_i)\sim y_i\f$,
 *  where Q is a quadratic polynomial.
 */
QuadraticPolynomialInterpolator(int nData, Vector x, Vector y);

/** Function value at t. */
double operator()(double t);

/** First derivative at t. */
double diff(double t);

/** Second derivative at t. */
double diff2(double t);


private:
	
	Vector a;  // coefficients Q(t) = a[0] + a[1]*t + a[2]*t*t 

};





#endif
