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
f1(Real t){ return 5*t*exp(-9*t*t/2); }

Real
f2(Real t){ return (1+t)*cos(7*PI*t); }

Real
f3(Real t){ return 5*t*exp(-9*t*t/2)+t*t*sin(11*PI*t); }

Real
f4(Real t){ return 1.6*abs(t)*sin(11*PI*t); }

Real
f5(Real t){ return t*t*(1+t); }


RealFunction
functionExample(int j)
{
   assert((0<=j)&&(j<=5));
   RealFunction f[6]={&f0,&f1,&f2,&f3,&f4,&f5};
   return f[j];
}



GPR_END_NAMESPACE(Gaussian)



