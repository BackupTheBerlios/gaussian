/***************************************************************************
                          RegressionPlots.h  -  description
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

#ifndef gpr_regressionplots_h
#define gpr_regressionplots_h

#include "gpr.h"


GPR_BEGIN_NAMESPACE(Gaussian)


/** Draws plots of regressions. */
class RegressionPlots {

public:

    /** The underlying Gaussian process regressor. */
    GPR& getGPR(){ return gpr_; }

    RegressionPlots(GPR& gpr): gpr_(gpr) { }


    /** PWD must contain directory "expansions".
     *  Plots regressions for all function examples 
     *  which are implemented. Regressors shown are \f$f_5,f_9,f_17\f$. 
     *  Plots are the directory "expansion/basis".
     */
    static void regressionPlots();


   /** Prints the values \f$\psi_i(t_j)\f$ of the current basis functions
    *  \f$\psi_i\f$ on [-1,+1] for \f$i=0,1,\dots,q\f$ at evenly spaced points
    *  \f$t_j\in[-1,+1]\f$, \f$j=0,1,...,m\f$ to the file BasisFunctionsData.txt in
	 *  gnuplot data format.
	 *
	 *  The first column contains the points \f$t_j\f$ and the next columns
	 *  contain \f$\psi_0(t_j),\psi_1(t_j),\dots,\psi_q(t_j)\f$. See
	 *  "doc/gnuplot_Readme.html" for instructions to plot these data
	 *  with gnuplot.
    **/
   void basisFunctionData(int q, int m);

private:

    GPR gpr_;
    
    // Sets the function underlying the data to example fj and writes plots
    // for fj with various values of N,n.
    void regressionPlots(int j, bool noisy);

};


GPR_END_NAMESPACE(Gaussian)


#endif

