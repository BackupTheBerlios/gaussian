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

};


/** Orthonormal Legendre polynomials on [-1+1] with normalized Lebesgue measure. **/
class LegendreBasis : public BasisFunctions {

public:

   RealArray1D values(Real x, int m);
   string name() const { return "Legendre_basis"; }
   Real ymin(int q){ return -5.0; }
   Real ymax(int q){ return 5.0; }
};


/** Orthonormal Fourier basis on [-1+1] with normalized Lebesgue measure. **/
class FourierBasis : public BasisFunctions {

public:

   RealArray1D values(Real x, int m);
   string name() const { return "Fourier_basis"; }
   Real ymin(int q){ return -2.2; }
   Real ymax(int q){ return 2.2; }
};
