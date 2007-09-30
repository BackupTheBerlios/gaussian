/***************************************************************************
 *   Copyright (C) 2007 by mjhmeyer   *
 *   mjhmeyer@yahoo.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#pragma once
#include "Typedefs.h"



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
CubicSpline(int nData, Vector x, Vector y, double ypLeft, double ypRight);
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