/***************************************************************************
                          Programs.cpp  -  description
                             -------------------
    begin                : Fri Jun 18 2004
    copyright            : (C) 2004 by 
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

#include "Programs.h"
#include "Array.h"
#include "gpr.h"
#include "RegressionPlots.h"
#include "plot.h"
#include <fstream>
using std::ofstream;



void
plotTest()
{
   // output file
   ofstream fout("test_plot.eps");
   PSPlot plot(fout,-1.0,+1.0,-1.2,+1.2);
   RealArray1D sine(401);
   RealArray1D expfn(401);

   for(int t=0;t<=400;t++){
      Real x=-1.0+1.0*t/200;
      sine[t]=sin(TWO_PI*x);
      expfn[t]=exp(-2*x*x);
   }

   plot.addFunction(expfn);
   plot.addFunction(sine);
}



void
interactiveRegression()
{
   GPR& gpr=GPR::setUp();
   int N=gpr.get_N();
   gpr.expansionData(N);
}



void
regressionPlots()
{
   RegressionPlots::regressionPlots();
}



