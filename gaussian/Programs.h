/***************************************************************************
                          Programs.h  -  description
                             -------------------
    begin                : Fri Jun 18 2004
    copyright            : (C) 2004 by Michael J.Meyer
    email                :  spyqqqdia@yahoo.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef gpr_programs_h
#define gpr_programs_h


#include "TypedefsMacros.h"
 
GPR_BEGIN_NAMESPACE(Gaussian)

void plotTest();
void choleskyTiming();

/** Writes the data files "ExpansionData.txt" and "FunctionData.txt"
 *  for plotting the regressors \f$f_0,f_1,\dots,f_N\f$ of the function \f$f(t)\f$
 *  with gnuplot. User chooses from several sample functions.
 *  The mean of the Gaussian prior P can be set to the origin, the vector of
 *  empirical coefficients or the vector (r,r,...,r), with r a real number.
 *  This last option is included to investigate the dependence of the
 *  algorithm on the mean of P.
 *
 *  The coefficients are computed from the data points \f$y_j=f(s_j)\f$,
 *  \f$0\leq j\leq n\f$. The points \f$s_j\f$ can be evenly spaced in [-1,+1]
 *  with \f$s_0=-1\f$, \f$s_n=+1\f$ or they can be uniformly random.
 *  The function values \f$y_j=f(s_j)\f$ can be exact or corrupted by
 *  independent Gaussian noise. Noisy data are of the form
 *  \f[y_j=f(s_j)+\sigma Z_j,\f]
 *  where the \f$Z_j\f$ are independent and standard normal. That is the standard
 *  deviation of the error is independent of the size of \f$f(s_j)\f$.
 *  The user supplies all parameters.
 **/
void interactiveRegression();


/** PWD must contain directory "expansions".
 *  Plots regressions for all function examples in all bases
 *  which are implemented. Regressors shown are \f$f_5,f_9,f_17\f$.
 *  Plots are the directory "expansion".
 */
void regressionPlots();


/** Estimates the integral and evaluation functional at t=0.5 for each of the
 *  function examples. The functionals are viewed both as linear functionals
 *  and general nonlinear functionals and the computation cartried out for both
 *  cases.
 */
void functionalEstimationTest();


/** Estimates the integral over [-1,+1] viewed as a linear functional
 *  and compares the value to the Monte Carlo
 *  integral using the data points and to the "exact" integral (Simpson's rule
 *  using 10000 equidistant points). Output written to file "Integrals.txt".
 *
 *  @param sigma standard deviation of data noise
 */
void bayesMonteCarlo(Real sigma);
 


GPR_END_NAMESPACE(Gaussian)

#endif
