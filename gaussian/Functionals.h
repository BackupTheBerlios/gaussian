/***************************************************************************
                          Functionals.h  -  description
                             -------------------
    begin                : Thu Jul 8 2004
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

#include <cassert>
#include "Array.h"

/*! \file Functionals.h
 *  Linear or nonlinear functionals
 *  \f[L\ :\ H=span{\ths\psi_0,\psi_1,\dots,\psi_N\ths\}\longrightarrow R\f]
 *  where H is the finite dimensional Hilbert space spanned by the basis
 *  used in a Gaussian process regressor {@link GPR}.
 *
 *  The class {@link LinearFunctional} derives from and implements all the 
 *  methods of {@link Functional} although linear functionals are treated
 *  completely differently in point estimation. We do this for cross validation
 *  of the methods.
 */


class GPR;


/** A general (not necessarily linear) functional $L$ on the span
 *  \f[H=span\{\ths\psi_0,\psi_1,\dots,\psi_N\ths\}\f]
 *  of the basis functions of a Gaussian process regressor.
 */
class Functional {

public:

   Functional(GPR* gp): gpr(gp) {  }

   /** The Gaussian process regressor defining the Hilbert space on
    *  which the functional is defined.
    */
   GPR* getGPR(){ return gpr; }

   /** The unconditional mean E(L), L=this, with respect to the measure
    *  underlying gpr.
    */
   virtual Real mean() const =0;

   /** The covariance \f$Cov(E_j,L)\f$, where L=this, \f$E_j=E_{s_j}\f$
    *  computed in the probability underlying gpr.
    */
   virtual Real covariance(int j) const =0;

private:

   GPR* gpr;

};


/** A linear functional $L$ on the span
 *  \f[H=span\{\ths\psi_0,\psi_1,\dots,\psi_N\ths\}\f]
 *  of the basis functions of a Gaussian process regressor.
 */
class LinearFunctional : public Functional {

public:

   LinearFunctional(GPR* gp): Functional(gp) {  }
   Real mean() const;
   Real covariance(int j) const;

   /** The sequence of values \f$L(\psi_k\f,\ k\leq N\f$, where the
    *  \f$\psi_k\f$ are the basis functions of gpr.
    */
   virtual const RealArray1D valuesOnBasisFunctions() const =0;


};


 
/** The integral \f$L(f)=\int_{-1}^{+1}f(t)dt\f$ as a linear as well as a general 
 *  functional. Used for "Bayesian Monte Carlo", see gpr-notes, section 6, p. 33.
 */
class Integral : public LinearFunctional {

public:

   Integral(GPR* gp): LinearFunctional(gp) {  }
   Real mean() const;
   Real covariance(int j) const;
   const RealArray1D valuesOnBasisFunctions() const;
};



/** The evaluation functional \f$L(f)=f(t)\f$ as a linear as well as a general
 *  functional. Used for point prediction and testing of the functional estimation
 *  methods.
 */
class EvaluationFunctional : public LinearFunctional {

public:

   /** Evaluation at the point s on the Hilbert space of gp. */
   EvaluationFunctional(Real s):
   LinearFunctional(gp);
   t(s)
   { assert((s<=1.0)&&(-1.0<=s)); }

   Real mean() const;
   Real covariance(int j) const;
   const RealArray1D valuesOnBasisFunctions() const;

private:

   Real t;   // point in [-1,+] where we evaluate

};


  