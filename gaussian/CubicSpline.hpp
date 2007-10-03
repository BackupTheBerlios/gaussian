#ifndef ____CubicSpline_H_161671___
#define ____CubicSpline_H_161671___

#pragma once
#include "Typedefs.hpp"



namespace spline {

/**!
*/
class CubicSpline {

public:

/** Cubic spline passing through the points (x[i],y[i]), 0<=i<nData.
 *  @param nData = length(x) = length(y) nuber of data points.
 *  @param ypLeft left endpoint derivative.
 *  @param ypRight right endpoint derivative.
 */
CubicSpline(int n_Data, Vector x, Vector y, double ypLeft, double ypRight);
~CubicSpline(void);
/** Number of knots. */
int nData(){ return nData_; }
/** Value at x*/
double operator()(double x){ return value(x);}
/** Value at x*/
double value(double x);
/** Derivative at x. */
double diff(double x);
/** Second derivative at x. */
double diff2(double x);
double x(int i) const { return x_[i]; }
double y(int i) const { return y_[i]; }
double h(int i) const { return h_[i]; }
/** Left endpoint derivative, read-write access. */
double& ypLeft(){ return ypLeft_; }
/** Right endpoint derivative, read-write access. */
double& ypRight(){ return ypRight_; }

private:
 
// i such that t\in(x_i,x_{i+1}]
int find(double t);
// read-write access
double& h(int i){ return h_[i]; }
double& g(int i){ return g_[i]; }
double& sp(int i){ return sp_[i]; }
  
		 int nData_;
		 Vector x_;
		 Vector h_;
		 Vector y_;
		 double ypLeft_;   // y'(x_0)
		 double ypRight_;  // y'(x_l), l=nData-1
		 Vector g_;         // gammas
		 Vector sp_;        // s'(x_{i+1})

};


} // end namespace spline


#endif
