/***************************************************************************
                          SampleFunctions.h  -  description
                             -------------------
    begin                : Tue Jun 15 2004
    copyright            : (C) 2004 by Michael J. Meyer
    email                : 
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef FunctionExamples_h
#define FunctionExamples_h

#include "TypedefsMacros.h"
#include "Array.h"


GPR_BEGIN_NAMESPACE(Gaussian)



/** \f$f(t)=sin(2\pi t)\f$.
 **/
Real f0(Real t);

/** \f$f(t)=5t*exp(-9t^2/2)\f$.
 **/
Real f1(Real t);

/** \f$f(t)=(1+t)cos(8\pi t)\f$.
 **/
Real f2(Real t);

/** \f$f(t)=f1(t)+t^2sin(11\pi t)\f$.
 **/
Real f3(Real t);

/** \f$f(t)=1.6|t|sin(10\pi t)\f$.
 **/
Real f4(Real t);

/** \f$f(t)=t^2+t^3\f$ */
Real f5(Real t);

/** The j-th function exmple. */
RealFunction functionExample(int j);


GPR_END_NAMESPACE(Gaussian)


#endif


