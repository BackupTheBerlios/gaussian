/***************************************************************************
                          Functionals.cpp  -  description
                             -------------------
    begin                : Thu Jul 8 2004
    copyright            : (C) 2004 by Michael J. Meyer
    email                : spyqqqdia@yahoo.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "Functionals.h"
#include "GPR.h"


//----------------LINEAR FUNCTIONAL----------------------------

Real
LinearFunctional::
mean()
{
   GPR* gpr=getGPR();
   int N=gpr->get_N();
   // the vector of unconditional expectations E[A_k]:
   const RealArray1D& mu=gpr->getPriorMean();
   // the sequence of values integral(gpr->psi_k)
   const RealArray1D l=valuesOnBasisFunctions((gpr);
   Real I=0.0;
   for(int k=0;k<=N;k++) I+=mu[k]*l[k];
   return I;
}


Real
LinearFunctional::
covariance(int j)
{
   // see gpr-notes, equation (39), p34.
   GPR* gpr=getGPR();
   int N=gpr->get_N();
   // the sequence of values integral(gpr->psi_k)
   const RealArray1D l=valuesOnBasisFunctions(gpr);
   // the matrix psi_k(s_j)
   const RealMatrix& gpr->get_psi();
   Real cv=0.0;
   for(int k=0;k<=N;k++) cv+=psi(k,j)*l[k];
   return cv;
}



//------------------INTEGRAL-------------------------------




const RealArray1D
Integral::
valuesOnBasisFunctions()
{
   GPR* gpr=getGPR();
   int N=gpr->get_N();
   return gpr->getBasisFunctions()->integrals(N);
}


//------------------EVALUATION FUNCTIONAL--------------------------





const RealArray1D
EvaluationFunctional::
valuesOnBasisFunctions()
{
   GPR* gpr=getGPR();
   int N=gpr->get_N();
   return gpr->getBasisFunctions()->values(t,N);
}

