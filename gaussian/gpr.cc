/***************************************************************************
 *            gpr.cc
 *
 *  Sat May 29 16:54:15 2004
 *  Copyright  2004  cpp
 *  cpp@linux
 ****************************************************************************/

#include "gpr.h"
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include "Random.h"
#include "plot.h"
using std::ofstream;
using std::system;
using std::setprecision;
using namespace Martingale;



//-------------------SET/GET----------------------



void
GPR::
setFunctionData(RealArray1D& w)
{
   y=w;
   computeEmpiricalCoefficients();
   if(have_gaussian) computeGaussianCoefficients();
}



void
GPR::
setFunction(RealFunction g)
{
   f=g;
   for(int j=0;j<=n;j++) y[j]=g(s[j]);
   computeEmpiricalCoefficients();
   if(have_gaussian) computeGaussianCoefficients();
}



void
GPR::
addNoise(Real sigma)
{
   for(int j=0;j<=n;j++) y[j]+=sigma*Random::sTN();
   computeEmpiricalCoefficients();
   if(have_gaussian) computeGaussianCoefficients();
}



void
GPR::
setRegressionType(RegressionType rt)
{
   if((rt==GAUSSIAN)&&!haveGaussian()) initGaussian();
   regrType=rt;
}



void
GPR::
setPriorMean(const RealArray1D& mean)
{
   if(regrType==GAUSSIAN)
   { mu=mean; computeGaussianCoefficients(); }
}



void
GPR::
setPriorMeanToOrigin()
{
    if(regrType==GAUSSIAN)
    { RealArray1D origin(n+1); mu=origin; computeGaussianCoefficients(); }
}



void
GPR::
setPriorMeanToEmpiricalCoefficients()
{
    if(regrType==GAUSSIAN)
    { mu=empCoeff; computeGaussianCoefficients(); }
}



void
GPR::
setBasisFunctions(BasisFunctions* bFcns)
{
    delete basis; basis=bFcns;
    initBasis();
    if(haveGaussian()) initGaussian();
}
   





//----------------CONSTRUCTORS------------------------------------------------


GPR::
GPR(int dmax, RealArray1D& t, RealArray1D& w, BasisFunctions* bFcns, RegressionType rt):
have_gaussian(false),
regrType(rt),
N(dmax),
n(t.getDimension()-1),
s(t),
y(w),
f(NULL),
basis(bFcns),
basis_name(basis->name()),
psi(N+1,n+1,0,0),
K(n+1), 
R(n+1),
mu(dmax+1),
a(dmax+1),
empCoeff(dmax+1)
{
   assert(bFcns!=NULL);
   initBasis();
   computeEmpiricalCoefficients();
   if(rt==GAUSSIAN){
      initGaussian();
      computeGaussianCoefficients();
      have_gaussian=true;
   }      
}




GPR::
GPR
(int dmax, RealArray1D& t, RealFunction g, Real sigma,
 BasisFunctions* bFcns, RegressionType rt):
have_gaussian(false),
regrType(rt),
N(dmax),
n(t.getDimension()-1),
s(t),
y(n+1),
f(g),
basis(bFcns),
basis_name(basis->name()),
psi(N+1,n+1,0,0),
K(n+1),
R(n+1),
mu(dmax+1),
a(dmax+1),
empCoeff(dmax+1)
{
   assert(bFcns!=NULL);
   for(int j=0;j<=n;j++) y[j]=f(s[j])+sigma*Random::sTN();
   initBasis();
   computeEmpiricalCoefficients();
   if(rt==GAUSSIAN){
      initGaussian();
      computeGaussianCoefficients();
      have_gaussian=true;
   }
}



GPR::
GPR(int M, int m, RealFunction g, BasisFunctions* bFcns, bool random):
have_gaussian(false),
regrType(EMPIRICAL),
N(M),
n(m),
s(m+1),
y(m+1),
f(g),
basis(bFcns),
basis_name(basis->name()),
psi(M+1,m+1,0,0),
K(m+1),
R(m+1),
mu(M+1),
a(M+1),
empCoeff(M+1)
{
   // abscissas s_j
   Real u=0.0; if(random) u=0.6;
   // ensure incresing s_j
   s[0]=-1.0;
   for(int j=1;j<n;j++) s[j]=-1.0+(j+u*Random::U01())*2.0/n;
   s[n]=1.0;

   for(int j=0;j<=n;j++) y[j]=f(s[j]); 

   initBasis();
   computeEmpiricalCoefficients();
}



//------------------INITIALIZATION--------------------------


void
GPR::
computeEmpiricalCoefficients()
{
      for(int k=0;k<=N;k++) empCoeff[k]=EC(k);
}



void
GPR::
initBasis()
{
   // intialize psi=(psi_k(s_j))_{0<=k<=N, 0<=j<=n}
	for(int j=0;j<=n;j++){

	   RealArray1D P=basisFunctionValues(s[j],N);
	   for(int k=0;k<=N;k++) psi(k,j)=P[k];
	}
}



void
GPR::
initGaussian()
{
   // intialize psi=(psi_k(s_j))_{0<=k<=N, 0<=j<=n}
	for(int j=0;j<=n;j++){

	   RealArray1D P=basisFunctionValues(s[j],N);
	   for(int k=0;k<=N;k++) psi(k,j)=P[k];
	}

	// allocate K=(K(s_i,s_j))
	for(int i=0;i<=n;i++)
	for(int j=i;j<=n;j++){

		Real sum=0.0;
		for(int k=0;k<=N;k++) sum+=psi(k,i)*psi(k,j);
		K(i,j)=sum;
	}
	LTRRealMatrix& rho=K.ltrRoot();
	R=rho;
	delete &rho;
	computeGaussianCoefficients();
}



void
GPR::
computeGaussianCoefficients()
{
   for(int k=0;k<=N;k++) a[k]=EA(k);
}



//-----------THE COEFFICIENTS-----------------------------


Real
GPR::
EC(int k)
{
   assert((0<=k)&&(k<=N));

	Real sum=0.0;
	for(int j=0;j<=n;j++) sum+=y[j]*psi(k,j);
	return sum/(n+1);
}



Real
GPR::
EA(int k)
{
  assert((0<=k)&&(k<=N));
	// compute the r[j]=R_{n+1,j}, 0<=j<=n,
	// equations 22,23 p11, gprs.ps
	RealArray1D r(n+1);
	r[0]=psi(k,0)/R(0,0);
	for(int j=1;j<=n;j++){

		Real sum=0.0;
		for(int i=0;i<j;i++)sum+=r[i]*R(j,i);

		r[j]=(psi(k,j)-sum)/R(j,j);
	}

  // the means mu_E[j]=E^P(E_j) of the evaluation functionals E_j
  // these are nonzero if P is not centered at zero
  RealArray1D mu_E(n+1);
  for(int j=0;j<=n;j++){

    Real sum=0.0;
    for(int k=0;k<=N;k++) sum+=mu[k]*psi(k,j);
    mu_E[j]=sum;
  }

	// compute the Z_j, 0<=j<=n
	RealArray1D Z(n+1);
	Z[0]=(y[0]-mu_E[0])/R(0,0);
	for(int j=1;j<=n;j++){

		Real sum=0.0;
		for(int k=0;k<j;k++) sum+=R(j,k)*Z[k];

		Z[j]=(y[j]-mu_E[j]-sum)/R(j,j);
	}

	// compute a_k=E(A_k), formula 20, p10, gprs.ps
  // note A_k has (unconditional) mean mu[k]
	Real a_k=mu[k];
	for(int j=0;j<=n;j++) a_k+=r[j]*Z[j];

	return a_k;
}



//-----------------EXPANSIONS-----------------------------------------


Real
GPR::
expansion(Real t, int q)
{
	RealArray1D P=basisFunctionValues(t,q);
	Real f_q=0.0;
   if(regrType==GAUSSIAN)
	   for(int k=0;k<=q;k++) f_q+=a[k]*P[k];
   else
      // empirical regression
      for(int k=0;k<=q;k++) f_q+=empCoeff[k]*P[k];
	return f_q;
}



void
GPR::
expansionData(int q)
{
   assert(q<=N);
   // write the data
	ofstream dout("FunctionData.txt");
	for(int j=0;j<=n;j++) dout << s[j] << "  " << y[j] << endl;
	dout.close();

	// write the expansions f_0,f_1,...,f_q evaluated at 801 points
	// to a file in gnuplot data format.
	ofstream fout("ExpansionData.txt");
	int m=800;
	for(int j=0;j<=m;j++){

	   Real t=-1.0+j*2.0/m;
	   fout << t << "  ";
      // if f is known (f!=NULL) write f
      if(f) fout << f(t) << "  ";
	   RealArray1D psi=basisFunctionValues(t,q);
	   Real sum=0.0;
      if(regrType==EMPIRICAL)
	     for(int i=0;i<=q;i++){ sum+=empCoeff[i]*psi[i]; fout << sum << "  "; }
      else
        for(int i=0;i<=q;i++){ sum+=a[i]*psi[i]; fout << sum << "  "; }
	   fout << endl;
   }
   fout.close();
}

			
//--------------------------TESTS------------------------------------		

void
GPR::
orthoTest(int q, int m)
{
    // write the matrix B=(P_i(s_j))_{0<=i<=q,0<=j<=m}
	// s_j\in[-1,+1] evenly spaced.
	RealMatrix B(q+1,m+1,0,0);
	for(int j=0;j<=m;j++){

	   Real x=-1.0+j*2.0/m;
	   RealArray1D P=basisFunctionValues(x,q);
	   for(int i=0;i<=q;i++) B(i,j)=P[i];
	}
	// Scale B by 1/sqrt(n), this scales BB' by 1/n and then
	// orthonormality with Monte Carlo integration is equivalent to BB'=I.
	Real f=1.0/sqrt(m+1); B*=f;

	cout << "Orthonormality: matrix of L^2 inner products (psi_i,psi_j):"
        << endl << endl
	     << B.aat();
}



void
GPR::
printBasisFunctions(int q)
{
   system("rm -R basis_plots");
   system("mkdir basis_plots");
   LegendreBasis().print(q);
   FourierBasis().print(q);
   system("mv *_basis basis_plots");
}




//-----------------USER SET UP---------------------------------

GPR&
GPR::
setUp()
{
   int N,n;
   RealArray1D t(1);
   Real sigma;
   RegressionType rType;
   BasisFunctions* bFcns;

   int ctr;   // prior mean centered? (0,1,2)
   Real r;    // 2: mu=(r,r,...,r)

   // DATA
   RealFunction  f;
   cout << "Sample function f:" << endl
        << "f0(t)=sin(2pi*t)...........[0]" << endl
        << "f1(t)=5t*exp(-9t^2/2)......[1]" << endl
        << "f2(t)=|t|..................[2]" << endl
        << "f3(t)=|t|^{1/3}............[3]" << endl << endl
        << "Enter f=[0,1,2,3]=";
   int fnum; cin>>fnum;
   // *f is reference to pointer to sample function
   switch(fnum){
      case 1  : f=&f1; break;
      case 2  : f=&f2; break;
      case 3  : f=&f3; break;
      default : f=&f0;
   }

   // noise level
   sigma=0.0;
   cout << "Noisy data (0/1)? Noisy=";
   int noisy; cin >> noisy;
   if(noisy==1){
      cout << "Standard deviation sigma=";
      cin>>sigma;
   }

   cout << endl << "Enter number N+1 of basis functions, N=";
	cin>>N;
   cout << "Enter number n+1 of data points, n=";
	cin>>n;
	t.resize(n+1);

	int random;
	cout << "Data points random (random=1) or evenly spaced (random=0), random=";
	cin>>random;

	// data abscissas s_j
	if(random==0)
	   for(int j=0;j<=n;j++) t[j]=-1.0+j*2.0/n;
   else
	   for(int j=0;j<=n;j++) t[j]=-1.0+2.0*Random::U01();


   // BASIS
   cout << "Which basis:" << endl
        << "Legendre basis..............[0]" << endl
        << "Fourier basis...............[1]" << endl
        << "Default is [0]" << endl << endl
        << "Basis=";
   int basis; cin >> basis;
   switch(basis){

      case 1  : bFcns=new FourierBasis(); break;
      default : bFcns=new LegendreBasis();
   }


   // REGRESSION TYPE
   cout << "Gaussian regression, enter parameters."
        << endl << endl
        << "Type of regression: " << endl
        << "Empirical................[0]" << endl
        << "Gaussian.................[1]" << endl
        << "Default is [0]" << endl << endl
        << "Regression=";
   int rt; cin>>rt;
   if(rt==1) rType=GPR::GAUSSIAN; else rType=GPR::EMPIRICAL;

   // REGRESSOR
   GPR* gpr=new GPR(N,t,f,sigma,bFcns,rType);

   // SET PRIOR MEAN
   if(rType==GPR::GAUSSIAN){

      // get center of Gaussian prior
      cout << "\nMean of Gaussian prior P:" << endl
           << "origin....................[0]" << endl
           << "empirical coefficients....[1]" << endl
           << "vector (r,r,...,r)........[2]"
           << "Default is [0]"
           << endl << endl
           << "Mean=[0,1,2]=";
      cin >> ctr;
      if(ctr==1) gpr->setPriorMeanToEmpiricalCoefficients();
      if(ctr==2){

         cout << "r="; cin>>r;
         RealArray1D nu(N+1);
         for(int k=0;k<=N;k++) nu[k]=r;
         gpr->setPriorMean(nu);
      }
   }

	// LOG FILE
	ofstream lout("ExpansionLog.txt");
	lout << setprecision(3);
	lout << "Function expansion, session log." << endl << endl
	     << n+1 << " data points, ";
	if(random==1)
		lout << "randomly spaced." << endl;
	else
		lout << "evenly spaced." << endl;

	if(noisy==1)
		lout << "Function data noisy, "
	        << "standard deviation of noise" << sigma << endl;
	else
		lout << "Function data exact." << endl;

   switch(fnum){
      case 1:  lout << "Function f(t)=5t*exp(-9t^2/2)" << endl; break;
      case 2:  lout << "Function f(t)=|t|" << endl; break;
      case 3:  lout << "Function f(t)=|t|^{1/3}" << endl; break;
      default: lout << "Function f(t)=sin(2pi*t)" << endl;
   }

   if(rType==GPR::GAUSSIAN){
		lout << "Type of regression: GAUSSIAN" << endl;
      if (ctr==1)
         lout << "Prior mean: empirical coefficients" << endl;
      else if (ctr==2)
         lout << "Prior mean: (r,r,...,r) with r=" << r;
      else
         lout << "Prior mean: origin" << endl;
   }
	else
		lout << "Type of regression: EMPIRICAL" << endl;

   switch(basis){

      case 1  : lout << "Basis: Fourier basis";
      default : lout << "Basis: Legendre polynomials";
   }

	lout.close();
   return *gpr;

} // setUp





