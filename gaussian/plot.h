/***************************************************************************
                          plot.h  -  description
                             -------------------
    begin                : Wed Jun 16 2004
    copyright            : (C) 2004 by Michael J. Meyer
    email                :  spyqqqdia2yahoo.com
    
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


 
#ifndef gaussian_plot_h
#define gaussian_plot_h

#include <plotter.h>
#include <fstream>
#include "Array.h"
using std::ofstream;
using namespace Martingale;


/** A postscript plotter. **/
class PSPlot : public PSPlotter {

   Real xmin, xmax, ymin, ymax;

   // plot cycles through 10 colors
   static const char* colors[10];
   int c;     // current color index

public:

   /** We must set the Plotter parameters
    *  <p><code>
    *  char* page="letter"; Plotter::parampl ("PAGESIZE",page);
    *  </code></p>
    *  before calling this constructor.
    **/
   PSPlot(ofstream& fout, Real x_min, Real x_max, Real y_min, Real y_max);

   ~PSPlot();

   /** The next pencolor in a cycle of 10. */
   const char* nextColor(){ c++; c=c%10; return colors[c]; }

   /** Adds the graph \f$y_j=f(t_j)\f$ with the \f$t_j\f$ evenly spaced
    *  in [xmin,xmax] and containing xmin, xmax.
    *
    * @param linewidth plotutils, plotter::flinewidth.
    * @param pencolor plotutils, plotter::pencolorname.
    **/
   void addFunction(Real* y, int n, Real linewidth, const char* pencolor);

   /** Adds the graph \f$y_j=f(t_j)\f$ with $n$ evenly spaced points \f$t_j\f$
    *  in [xmin,xmax] containing xmin, xmax.
    **/
   void addFunction(Real* y, int n);

   /** Adds the graph \f$y_j=f(t_j)\f$ with the \f$t_j\f$ evenly spaced
    *  in [xmin,xmax] and containing xmin, xmax.
    *
    * @param linewidth plotutils, plotter::flinewidth.
    * @param pencolor plotutils, plotter::pencolorname.
    **/
   void addFunction(RealArray1D& y, Real linewidth, const char* pencolor);

   /** Adds the graph \f$y_j=f(t_j)\f$ with the \f$t_j\f$ evenly spaced
    *  in [xmin,xmax] and containing xmin, xmax.
    **/
   void addFunction(RealArray1D& y);

   /** Adds the points \f$(t_j,y_j)\f$ as solid red squares.
    *  y and t must be arrays of at least length n.
    **/
   void addPoints(Real* y, Real* t, int n);

   /** Adds the points \f$(t_j,y_j)\f$ as solid red squares.
    *  y and t must be arrays of the same length.
    **/
   void addPoints(const RealArray1D& y, const RealArray1D& t);

   /** Draw text label at (x,y) */
   void drawLabel(Real x, Real y, string s);

   /** Draws diamond shaped marker centered at (x,y). */
   void drawMarker(Real x, Real y);

}; 



#endif


