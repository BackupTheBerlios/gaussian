
#include "Utils.h"
#include <plot.h>
#include <fstream>
#include <cmath>
#include <cstdlib>
using std::ofstream;
using std::sin;
using std::exp;
using std::min;
using std::max;



const char*
PSPlot::
colors[10]={ "brown", "black", "red", "green", "blue", "orange", 
              "yellow", "magenta", "cyan", "darkblue" };



PSPlot::
PSPlot(ofstream& fout, Real x_min, Real x_max, Real y_min, Real y_max):
PSPlotter(cin,fout,cerr),
xmin(x_min),
xmax(x_max),
ymin(y_min),
ymax(y_max),
c(0)
{
    if (openpl () < 0)            
    {  cerr << "Couldn't open Plotter\n"; exit(EXIT_FAILURE); }

    // margins around bounding box
    Real xm=(xmax-xmin)/15,
         ym=(ymax-ymin)/15,
         tmin=xmin-xm, tmax=xmax+xm,
         wmin=ymin-ym, wmax=ymax+ym,
         wspan=wmax-wmin;

    // display is square, so set fspace so that bounding box has
    // golden ratio
    fspace(tmin,wmin,tmax,wmin+(1+GR)*wspan);    
     
    // draw biggest box displayed
    flinewidth(0.0);
    fbox(tmin,wmin,tmax,wmax);     
    pencolorname("black"); 
    fbox(xmin,ymin,xmax,ymax);      
}


PSPlot::
~PSPlot()
{
    if (closepl () < 0)
    {  cerr << "Couldn't close Plotter\n"; }
}


void
PSPlot::
addFunction(Real* y, int n, Real lw, const char* pc)
{
    flinewidth(lw);
    pencolorname(pc);

    Real x=xmin, dx=(xmax-xmin)/(n-1);
    x=xmin; fmove(x,y[0]);
    for(int i=1;i<n;i++){ x+=dx; fcont(x,y[i]); }
    endpath();
}


void
PSPlot::
addFunction(Real* y, int n)
{
   addFunction(y,n,0.0,nextColor());
}


void
PSPlot::
addFunction(RealArray1D& y, Real lw, const char* pc)
{
    int n=y.getDimension();
    addFunction(y.getData(),n,lw,pc);
}

   
void
PSPlot::
addFunction(RealArray1D& y)
{
    addFunction(y,0.0,nextColor());
}









   
      