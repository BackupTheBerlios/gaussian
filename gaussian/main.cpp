/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Sun Jun  6 03:04:25 PDT 2004
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Programs.h"
#include "gpr.h"
#include <cstdlib>

using namespace Gaussian;

int main(int argc, char *argv[])
{
  //interactiveRegression();
  //GPR::printBasisFunctions(9);
  //regressionPlots();
  //plotTest();
  //choleskyTiming();
  functionalEstimationTest();

  return EXIT_SUCCESS;
}
