/***************************************************************************
                          RegressionPlots.cpp  -  description
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



#include "RegressionPlots.h"
#include <fstream>
#include <iomanip>       // setprecision
#include <cstdlib>
#include <cassert>
#include "Random.h"
#include "plot.h"
#include "Utils.h"
using std::ofstream;
using std::system;
using namespace Martingale;








void
RegressionPlots::
basisFunctionData(int q, int m)
{
	ofstream fout("BasisFunctionData.txt");

	for(int j=0;j<=m;j++){

      Real t=-1.0+j*2.0/m;
	   fout << t << "  ";
      RealArray1D psi_t=gpr_.basisFunctionValues(t,q);
	   for(int k=0;k<=q;k++) fout << psi_t[k] << "  ";
	   fout << endl;
   }
   fout.close();
}



void
RegressionPlots::
regressionPlots(int j, bool noisy)
{
    int n=gpr_.get_n();
    int N=gpr_.get_N();
    const RealArray1D& y=gpr_.getFunctionData();
    const RealArray1D& s=gpr_.get_s();
    
    RealFunction f=functionExample(j);   
    gpr_.setFunction(f);
    if(noisy) gpr_.addNoise(0.2);

    // allocate the directories)
    string basis=gpr_.getBasisFunctions()->name();
    // for string conversion
    Int_t j_t(j);
    Int_t n_t(n);
    // allocate directories
    // fj/Ecact
    string fcn="f"+j_t.str();
    string fcn_dir="expansions/"+basis+"/"+fcn;

    if(!noisy){

       system(("mkdir "+fcn_dir).c_str());
       // fj/Exact
       system(("mkdir "+fcn_dir+"/Exact").c_str());
       system(("mkdir "+fcn_dir+"/Exact/Empirical").c_str());
       system(("mkdir "+fcn_dir+"/Exact/Gaussian").c_str());
    }else{
       system(("mkdir "+fcn_dir+"/Noisy").c_str());
       system(("mkdir "+fcn_dir+"/Noisy/Empirical").c_str());
       system(("mkdir "+fcn_dir+"/Noisy/Gaussian").c_str());
    }

    // expansions
    RealArray1D f_true(801);
    RealArray1D f_q(801);

    // EMPIRICAL EXPANSIONS
    gpr_.setRegressionType(GPR::EMPIRICAL);
    int q=gpr_.polluteAndPredictCV();
    Int_t q_t(q);
    string q_label="q = "+q_t.str();
    // plotter
    string file_e=fcn+"-"+n_t.str();
    if(noisy) file_e+="-0.2"; file_e+="-emp.eps";
    ofstream fout_e(file_e.c_str());
    PSPlot plot_e(fout_e,-1.0,1.0,-2.0,2.0);


    for(int j=0;j<=800;j++){

       Real t=-1.0+j*1.0/400;
       f_true[j]=f(t);
       f_q[j]=gpr_.expansion(t,q);
    }
    plot_e.addFunction(f_true);        // green
    plot_e.addFunction(f_q);           // blue
    plot_e.addPoints(y,s);
    plot_e.drawLabel(-0.9,-1.9,q_label);

    // GAUSSIAN EXPANSIONS
    gpr_.setRegressionType(GPR::GAUSSIAN);
    q=gpr_.polluteAndPredictCV();
    Int_t qt(q);
    q_label="q = "+qt.str();
    // plotter
    string file_g=fcn+"-"+n_t.str();
    if(noisy) file_g+="-0.2"; file_g+="-gsn.eps";
    ofstream fout_g(file_g.c_str());
    PSPlot plot_g(fout_g,-1.0,1.0,-2.0,2.0);


    for(int j=0;j<=800;j++){

       Real t=-1.0+j*1.0/400;
       f_true[j]=f(t);
       f_q[j]=gpr_.expansion(t,q);
    }
    plot_g.addFunction(f_true);    
    plot_g.addFunction(f_q);
    plot_g.addPoints(y,s);
    plot_g.drawLabel(-0.9,-1.9,q_label);

    // move the files into appropriate directories
    if(noisy){
       
       string cmd_e="mv "+file_e+" expansions/"+basis+"/"+fcn+"/Noisy/Empirical/";
       string cmd_g="mv "+file_g+" expansions/"+basis+"/"+fcn+"/Noisy/Gaussian/";
       system(cmd_e.c_str());
       system(cmd_g.c_str());
    } else {
       string cmd_e="mv "+file_e+" expansions/"+basis+"/"+fcn+"/Exact/Empirical/";
       string cmd_g="mv "+file_g+" expansions/"+basis+"/"+fcn+"/Exact/Gaussian/";
       system(cmd_e.c_str());
       system(cmd_g.c_str());
    }       

} // end regressionPlots




void
RegressionPlots::
regressionPlots()
{
   // BASIS
   BasisFunctions* bFcns;
   cout << "Which basis:" << endl
        << "Legendre basis..............[0]" << endl
        << "Fourier basis...............[1]" << endl
        << "Default is [0]" << endl << endl
        << "Basis=";
   int bs; cin >> bs;
   switch(bs){

      case 1  : bFcns=new FourierBasis(); break;
      default : bFcns=new LegendreBasis();
   }
   // design model
   cout << "Data points s_j evenly spaced (0) or random (1); spacing=";
   int spcng; cin>>spcng;
   bool random=false; if(spcng==1) random=true;

   string basis=bFcns->name();
   string cmd="mkdir expansions/"+basis;
   system(cmd.c_str());

   for(int j=0;j<5;j++)
   for(int n=10;n<91;n+=10){

      int N=n; // number of basis functions
      GPR gpr(N,n,functionExample(j),bFcns,random); // clean data
      RegressionPlots regrPlot(gpr);
      bool noisy=false;
      regrPlot.regressionPlots(j,noisy);
      regrPlot.getGPR().addNoise(0.2);
      noisy=true;      
      regrPlot.regressionPlots(j,noisy);
   }
}





