/***************************************************************************
                          BasisFunctions.h  -  description
                             -------------------
    begin                : Tue Jun 15 2004
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


 #include "Matrix.h"
 #include <string>
 using std::string;
 using namespace Martingale;

 
/** Basis functions on [-1,+1]. Usually orthonormal in \f$L^2([-1,+1], dt)\f$.
 *  Note that we are using <i>normalized</i> Lebesgue measure \f$d\mu=(1/2)dt\f$
 *  on [-1+1].
 **/
class BasisFunctions {

public:

   /** Computes the sequence \f$(\psi_0(x),\dots,\psi_m(x))\f$ of values of the
    *  basis functions.
    **/
   virtual RealArray1D values(Real x, int m)=0;

   /** Multiplier \f$\rho\in(0,1]\f$ used in {@link GPR::polluteAndPredictCV()}.
    *  We move from \f$f_{q-1}\f$ to \f$f_q\f$ only if error drops of by factor
    *  \f$\rho\f$.
    **/
   virtual Real roughnessPenalty(int q){ return 1.0; }
   
   /** Designation of basis.
    **/
   virtual string name() const =0;

   /** Designation of basis.
    **/
  const char* cname() const { return name().c_str(); }

   /** Tight lower bound for the first q basis functions.
    *  Used only in drawing graphs of the basis functions.
    */
   virtual Real ymin(int q)=0;

   /** Tight upper bound for the first q basis functions.
    *  Used only in drawing graphs of the basis functions.
    */
   virtual Real ymax(int q)=0;   

   /** Print an EPS graph of the first q basis functions. **/
   void print(int q);

   /** The sequence of integrals
    *  \f$I_k=\int_{-1}^{+1}\psi_k(t)dt,\ k\leq N.\f$
    */
   virtual const RealArray1D integrals(int N);

};


/** Orthonormal Legendre polynomials on [-1+1] with normalized Lebesgue measure. **/
class LegendreBasis : public BasisFunctions {

public:

   RealArray1D values(Real x, int m);
   Real roughnessPenalty(int q);
   string name() const { return "Legendre_basis"; }
   Real ymin(int q){ return -5.0; }
   Real ymax(int q){ return 5.0; }
   const RealArray1D integrals(int N);
};


/** Orthonormal Fourier basis on [-1+1] with normalized Lebesgue measure. **/
class FourierBasis : public BasisFunctions {

public:

   RealArray1D values(Real x, int m);
   Real roughnessPenalty(int q);
   string name() const { return "Fourier_basis"; }
   Real ymin(int q){ return -2.2; }
   Real ymax(int q){ return 2.2; }
   const RealArray1D integrals(int N);
};
