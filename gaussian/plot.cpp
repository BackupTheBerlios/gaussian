
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


GPR_BEGIN_NAMESPACE(Gaussian)




const char*
PSPlot::
colors[10]={ "brown", "green", "blue", "orange", "red", "black",
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
    // golden ratio (horizontal squeeze by 1+GR).
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


void
PSPlot::
addPoints(Real* y, Real* t, int n)
{
   for(int j=0;j<n;j++) drawMarker(t[j],y[j]);
}


void
PSPlot::
addPoints(const RealArray1D& y, const RealArray1D& t)
{
   int n=y.getDimension(),
       m=t.getDimension();
   assert(n==m);
   Real* w=y.getData();
   Real* s=t.getData();
   addPoints(w,s,n);
}


void
PSPlot::
drawLabel(Real x, Real y, string s)
{
    Real dx=xmax-xmin,
         dy=ymax-ymin;
    ffontsize(0.07);
    pencolorname("black");
    fscale(1.0,(1.0+GR)*dy/dx);
    fmove(x,y*dx/((1.0+GR)*dy));
    // undo horizontal distortion
    label(s.c_str());
    fscale(1.0,dx/((1.0+GR)*dy));
}


void
PSPlot::
drawMarker(Real x, Real y)
{
   colorname("red");
   filltype(1);
   // compensate for the squeeze in dy by 1+GR   
   Real dx=(xmax-xmin)/400, dy=dx*(1+GR);

   // left corner of diamond
   fmove(x-dx,y);
   flinerel(0,0,dx,dy);
   flinerel(0,0,dx,-dy);
   flinerel(0,0,-dx,-dy);
   flinerel(0,0,-dx,dy);
   endpath();

   colorname("white");
   filltype(0);
}   



GPR_END_NAMESPACE(Gaussian)

   
      