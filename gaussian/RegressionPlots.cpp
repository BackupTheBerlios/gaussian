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
    
    RealFunction f=functionExample(j);   
    gpr_.setFunction(f);
    if(noisy) gpr_.addNoise(0.2);

    // allocate the directories)
    string basis=gpr_.getBasisFunctions()->name();
    // for string conversion
    Int_t j_t(j);
    Int_t n_t(n);
    Int_t N_t(N);
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
    RealArray1D f_5(801);
    RealArray1D f_9(801);
    RealArray1D f_17(801);

    // EMPIRICAL EXPANSIONS
    gpr_.setRegressionType(GPR::EMPIRICAL);
    // plotter
    string file_e=fcn+"-"+N_t.str()+"-"+n_t.str();
    if(noisy) file_e+="-0.2"; file_e+=+"-emp.eps";
    ofstream fout_e(file_e.c_str());
    PSPlot plot_e(fout_e,-1.0,1.0,-2.0,2.0);


    for(int j=0;j<=800;j++){

       Real t=-1.0+j*1.0/400;
       f_true[j]=f(t);
       f_5[j]=gpr_.expansion(t,5);
       f_9[j]=gpr_.expansion(t,9);
       f_17[j]=gpr_.expansion(t,17);
    }
    plot_e.addFunction(f_true);
    plot_e.addFunction(f_5);
    plot_e.addFunction(f_9);
    plot_e.addFunction(f_17);
cerr << endl << endl << endl << "f: " << f_true << endl << endl;
cerr << endl << endl << endl << "f_5: " << f_5 << endl << endl;
cerr << endl << endl << endl << "f_9: " << f_9 << endl << endl;

    // GAUSSIAN EXPANSIONS
    gpr_.setRegressionType(GPR::GAUSSIAN);
    // plotter
    string file_g=fcn+"-"+N_t.str()+"-"+n_t.str();
    if(noisy) file_g+="-0.2"; file_g+=+"-gsn.eps";
    ofstream fout_g(file_g.c_str());
    PSPlot plot_g(fout_g,-1.0,1.0,-2.0,2.0);


    for(int j=0;j<=800;j++){

       Real t=-1.0+j*1.0/400;
       f_true[j]=f(t);
       f_5[j]=gpr_.expansion(t,5);
       f_9[j]=gpr_.expansion(t,9);
       f_17[j]=gpr_.expansion(t,17);
    }
    plot_g.addFunction(f_true);
    plot_g.addFunction(f_5);
    plot_g.addFunction(f_9);
    plot_g.addFunction(f_17);

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

   for(int j=0;j<4;j++)
   for(int n=11;n<32;n+=10){

      int N=int(1.6*n);
      GPR gpr(N,n,functionExample(j),bFcns,random); // clean data
      RegressionPlots regrPlot(gpr);
      bool noisy=false;
      regrPlot.regressionPlots(j,noisy);
      regrPlot.getGPR().addNoise(0.2);
      noisy=true;      
      regrPlot.regressionPlots(j,noisy);
   }
}





