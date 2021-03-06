/***************************************************************************
                          BasisFunctions.cpp  -  description
                             -------------------
    begin                : Wed Jun 16 2004
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


#include "BasisFunctions.h"
#include <cmath>
#include "plot.h"
using std::sin;
using std::cos;

GPR_BEGIN_NAMESPACE(Gaussian)


void
BasisFunctions::
print(int q)
{
   ofstream fout(cname());
   PSPlot plot(fout,-1.0,+1.0,ymin(q),ymax(q));
   // evaluate each basis function at 800 evenly spaced points in [-1,+1]
   // and fill the first q basis functions into the rows of psi
   RealArray2D psi(q,801);   
   for(int j=0;j<801;j++){

      Real t=-1.0+j*0.0025;
      RealArray1D psi_t=values(t,q);
      for(int k=0;k<q;k++) psi(k,j)=psi_t[k];
   }
   // add the rows to the plot
   for(int k=0;k<q;k++) plot.addFunction(psi[k],801);
}



RealArray1D
BasisFunctions::
integrals(int N) const
{
   cerr << "BasisFunctions: integrals not implemented current basis. Exiting.";
   exit(EXIT_FAILURE);
   RealArray1D I(N+1);
   return I;
}
      
      
      
//---------------- LEGENDRE BASIS ----------------------------      


RealArray1D
LegendreBasis::
values(Real x, int m) const
{
	RealArray1D P(m+1);
	int r;
	P[0]=1.0;
	P[1]=sqrt(3)*x;
	for(int n=1;n<m;n++){

		r=2*n+1;
		P[n+1]=sqrt(r*(r+2))*x*P[n];
		P[n+1]-=n*sqrt(1.0+4.0/(r-2))*P[n-1];
		P[n+1]/=(n+1);
	}
	return P;
}


Real
LegendreBasis::
roughnessPenalty(int q) const
{
   return exp(-(q/5)*0.3);
}


RealArray1D
LegendreBasis::
integrals(int N) const
{
   RealArray1D I(N+1);
   I[0]=2.0;    // all other integrals are zero
   return I;
}




//---------------- FOURIER BASIS ----------------------------  


RealArray1D
FourierBasis::
values(Real x, int m)  const
{
   assert(m>=0);
   RealArray1D P(m+1);
   // recall that we are using normalized Lebesgue measure 0.5*dt on [-1,+1]
   // (affects the empirical coefficients)
   for(int k=0;k<=m;k+=2) P[k]=SQRT_TWO*cos((k/2)*PI*x);
   for(int k=1;k<=m;k+=2) P[k]=SQRT_TWO*sin((1+k/2)*PI*x);
   
	return P;
}


Real
FourierBasis::
roughnessPenalty(int q) const
{
   return exp(-(q/5)*0.3);
}


RealArray1D
FourierBasis::
integrals(int N) const
{
   RealArray1D I(N+1);
   I[0]=2.0*SQRT_TWO;    // all other integrals are zero
   return I;
}


GPR_END_NAMESPACE(Gaussian)

