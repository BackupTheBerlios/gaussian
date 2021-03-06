/***************************************************************************
                          Programs.cpp  -  description
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

#include "Programs.h"
#include "Utils.h"
#include "Matrix.h"
#include "tnt/tnt/tnt_array2d.h"
#include "tnt/jama/jama_cholesky.h"
#include "gpr.h"
#include "RegressionPlots.h"
#include "plot.h"
#include "Functionals.h"
#include <fstream>
using std::ofstream;


GPR_BEGIN_NAMESPACE(Gaussian)


void
plotTest()
{
   // output file
   ofstream fout("test_plot.eps");
   PSPlot plot(fout,-1.0,+1.0,-2.2,+2.2);
   RealArray1D sine(401);
   RealArray1D expfn(401);
   RealArray1D t(21);
   RealArray1D y(21);

   for(int t=0;t<=400;t++){
      Real x=-1.0+1.0*t/200;
      sine[t]=f3(x);
      expfn[t]=exp(-2*x*x);
   }
   for(int j=0;j<=20;j++){

       t[j]=-1.0+j*0.1;
       y[j]=expfn[20*j];
   }

   plot.addFunction(expfn);
   plot.addFunction(sine);
   plot.addPoints(y,t);

   //label
   string label="q = 11";
   plot.drawLabel(-0.95,-1.95,label);
   cout << endl << endl
        << "Plot in file test_plot.eps.";
}


void choleskyTiming()
{
    bool done=false;
    while(!done){
    cout << "Dimension n=";
    int n; cin>>n;

    cout << endl << "Doing Martingale::Cholesky...";
    Timer watch;

    UTRRealMatrix K(n);
    for(int i=0;i<n;i++)
    for(int j=i;j<n;j++) K(i,j)=1.0*(1+i)/(1+j);

    watch.start();
    LTRRealMatrix R=K.ltrRoot();
    watch.stop();
    watch.report("Martingale::Cholesky:");

    TNT::Array2D<Real> M(n,n);
    TNT::Array2D<Real> L(n,n);
    for(int i=0;i<n;i++)
    for(int j=i;j<n;j++) M[i][j]=M[j][i]=1.0*(1+i)/(1+j);

    watch.start();
    JAMA::Cholesky<Real> chol(M);
    L=chol.getL();
    watch.stop();
    watch.report("TNT::Cholesky:");

    cout << endl << "Again (0/1)? Again=";
    cin>>n; if(n==0) done=true;
    } // end main loop
}    
    



void
interactiveRegression()
{
   bool done=false;
   while(!done){

     GPR& gpr=GPR::setUp();
     gpr.conditioning();
     int N=gpr.get_N();
     gpr.expansionData(N);

     cout << endl << endl << endl
          << "Do again? (0/1), again=";
     int again=0; cin>>again;
     if(again==0) done=true;
   }
}



void
regressionPlots()
{
   RegressionPlots::regressionPlots();
}


void
functionalEstimationTest()
{
   bool done=false;
   while(!done){

     GPR* gpr=&(GPR::setUp());
     gpr->conditioning();

     cout << endl << "Integral estimation (Bayesian Monte Carlo):";
     Integral I(gpr);
     Real I0=gpr->estimateFunctional(I);
     cout << endl << "Value as general functional = " << I0;

     I0=gpr->estimateLinearFunctional(I);
     cout << endl << "Value as linear functional = " << I0;

     I0=gpr->monteCarloIntegral();
     cout << endl << "Monte Carlo integral = " << I0;

     cout << endl << endl << "Point estimation at t=0.5:";
     EvaluationFunctional L(gpr,0.5);
     I0=gpr->estimateFunctional(L);
     cout << endl << "Value as general functional = " << I0;

     I0=gpr->estimateLinearFunctional(L);
     cout << endl << "Value as linear functional = " << I0;
     

     cout << endl << endl << endl
          << "Do again? (0/1), again=";
     int again=0; cin>>again;
     if(again==0) done=true;
   } // end while
}



void
writeIntegrals(ofstream& os, GPR* gpr)
{
   os << gpr->getBasisFunctions()->name() << endl << endl
      << "Function     Monte Carlo    Bayesian Monte Carlo    Exact integral"
      << endl << endl;
   for(int l=0;l<N_FUNC;l++){

      RealFunction f=functionExample(l);
      gpr->setFunction(f);
      Integral I(gpr);
      Real Simpson=integral(-1.0,1.0,f),
           MC=gpr->monteCarloIntegral(),
           BMC=gpr->estimateLinearFunctional(I);

      os << "  f" << l << "      " << MC << "        " << BMC
         << "          " << Simpson << endl;
   }
   os << endl << endl;
}



void
bayesMonteCarlo(Real sigma)
{
   cout << "Computing Bayesian Monte Carlo intergrals" << endl
        << "for all function examples in all bases:" << endl
        << "how many data points? n=";
   int n; cin >> n;
      // data points t
   int N=n;
   RealArray1D t=GPR::dataPoints(n,false);
   ofstream fout("Integrals.txt");
   fout << "Data points:" << endl << t << endl << endl;
   
   BasisFunctions* bFcns=new LegendreBasis();
   RealFunction f=functionExample(0);
   GPR* gpr=new GPR(N,t,f,bFcns,sigma,GPR::GAUSSIAN);
   gpr->conditioning();

   writeIntegrals(fout,gpr);

   bFcns=new FourierBasis();
   gpr->setBasisFunctions(bFcns);

   writeIntegrals(fout,gpr);

   fout.close();

   cout << endl << "Done, data in file Integrals.txt.";
}




GPR_END_NAMESPACE(Gaussian)
