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


#include "FunctionExamples.h"
#include "Random.h"
#include <cmath>
#include <fstream>
#include <iomanip>    // setprecision
#include <cassert>

using std::ofstream;
using std::setprecision;
using std::sin;
using std::exp;
using std::log;
using std::abs;


GPR_BEGIN_NAMESPACE(Gaussian)




Real
f0(Real t){ return sin(TWO_PI*t); }

Real
f1(Real t){ return 4*(1-t)*exp(-9*t*t/2); }

Real
f2(Real t){ return (1+t)/(1+t*t); }

Real
f3(Real t){ return (1+t+t*t)/(2-t+t*t); }

Real
f4(Real t){ return 1.6*abs(t)*sin(11*PI*t); }

Real
f5(Real t){ return t*t+exp(-t/2); }


RealFunction
functionExample(int j)
{
   assert((0<=j)&&(j<N_FUNC));
   RealFunction f[6]={&f0,&f1,&f2,&f3,&f4,&f5};
   return f[j];
}


Real
integral(Real a, Real b, RealFunction f)
{
   int N=10000;
   Real h=(b-a)/N;
   Real I0=f(a)+f(b), I1=0.0, I2=0.0;

   for(int i=1;i<N;i+=2){ Real t=a+i*h; I1+=f(t); }
   for(int i=2;i<N;i+=2){ Real t=a+i*h; I2+=f(t); }

   return h*(I0+2*I2+4*I1)/3;
} 




GPR_END_NAMESPACE(Gaussian)



