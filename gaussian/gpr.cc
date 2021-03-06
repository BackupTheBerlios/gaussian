/***************************************************************************
 *            gpr.cc
 *
 *  Sat May 29 16:54:15 2004
 *  Copyright  2004  cpp
 *  cpp@linux
 ****************************************************************************/

#include "gpr.h"
#include "Functionals.h"
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include "Random.h"
#include "plot.h"
#include <cmath>
#include <vector>
#include <algorithm>
using std::ofstream;
using std::system;
using std::setprecision;
using std::sqrt;
using std::vector;
using std::sort;


GPR_BEGIN_NAMESPACE(Gaussian)



Real
GPR::
EPS=0.000000001;



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
   sigma_=sqrt(sigma_*sigma_+sigma*sigma);
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
    { mu=b; computeGaussianCoefficients(); }
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
GPR(int dmax, RealArray1D& t, RealArray1D& w, BasisFunctions* bFcns,
    Real sigma, RegressionType rt):
have_gaussian(false),
regrType(rt),
N(dmax),
n(t.getDimension()-1),
sigma_(sigma),
s(t),
y(w),
yp(t.getDimension()),
f(NULL),
basis(bFcns),
basis_name(basis->name()),
psi(N+1,n+1,0,0),
K(n+1), 
R(n+1),
mu(dmax+1),
a(dmax+1),
b(dmax+1),
ap(dmax+1),
bp(dmax+1)
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
(int dmax, RealArray1D& t, RealFunction g, 
 BasisFunctions* bFcns, Real sigma, RegressionType rt):
have_gaussian(false),
regrType(rt),
N(dmax),
n(t.getDimension()-1),
sigma_(sigma),
s(t),
y(n+1),
yp(n+1),
f(g),
basis(bFcns),
basis_name(basis->name()),
psi(N+1,n+1,0,0),
K(n+1),
R(n+1),
mu(dmax+1),
a(dmax+1),
b(dmax+1),
ap(dmax+1),
bp(dmax+1)
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
sigma_(0.0),
s(m+1),
y(m+1),
yp(m+1),
f(g),
basis(bFcns),
basis_name(basis->name()),
psi(M+1,m+1,0,0),
K(m+1),
R(m+1),
mu(M+1),
a(M+1),
b(M+1),
ap(M+1),
bp(M+1)
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



//------------------FUNCTIONAL ESTIMATION--------------------------


Real
GPR::
estimateFunctional(Functional& L)
{
   // compute the r[j]=R_{n+1,j}, 0<=j<=n,
   // gpr-notes (36),(37), p. 34.
   RealArray1D r(n+1);
   r[0]=L.covariance(0)/R(0,0);
   for(int j=1;j<=n;j++){

      Real sum=0.0;
      for(int i=0;i<j;i++)sum+=r[i]*R(j,i);

      r[j]=(L.covariance(j)-sum)/R(j,j);
   }

   // the means mu_E[j]=E^P(E_j) of the evaluation functionals E_j
   // these are nonzero if P is not centered at zero
   RealArray1D mu_E(n+1);
   for(int j=0;j<=n;j++){

      Real sum=0.0;
      for(int k=0;k<=N;k++) sum+=mu[k]*psi(k,j);
      mu_E[j]=sum;
   }

   // compute the Z_j, 0<=j<=n, gpr-notes (25), p. 11 with the
   // (unconditional) means E(E_j) added on the right hand side.
   RealArray1D Z(n+1);
   Z[0]=(y[0]-mu_E[0])/R(0,0);
   for(int j=1;j<=n;j++){

      Real sum=0.0;
      for(int k=0;k<j;k++) sum+=R(j,k)*Z[k];

      Z[j]=(y[j]-mu_E[j]-sum)/R(j,j);
   }

   // compute EL=E[L|data], gpr-notes (25), p. 11 with the (unconditional) means
   // added on the right hand side.
   Real EL=L.mean();
   for(int j=0;j<=n;j++) EL+=r[j]*Z[j];

   return EL;
}



Real
GPR::
estimateLinearFunctional(LinearFunctional& L)
{
   RealArray1D d=getCoefficients();
   // the sequence of values L(psi_k), k<=N.
   RealArray1D l=L.valuesOnBasisFunctions();
   // use regressor f_N to preserve equality with general
   // functional estimation, gpr-notes, equation (40), p34.
   // Recall: L(this,N)[k]=L(psi_k), k<=N.
   Real EL=0.0;
   for(int k=0;k<=N;k++) EL+=d[k]*l[k];
   return EL;
}



Real
GPR::
monteCarloIntegral()
{
   Real I=0.0;
   for(int j=0;j<=n;j++) I+=y[j];
   return 2*I/(n+1);
}


//------------------INITIALIZATION--------------------------


void
GPR::
computeEmpiricalCoefficients()
{
      for(int k=0;k<=N;k++) b[k]=EC(k,y);
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
  // add SIGMA along the diagonal to make sure
  // this matrix is positive definite
  for(int i=0;i<=n;i++) K(i,i)+=SIGMA;
  
  // Cholesky root R
	LTRRealMatrix& rho=K.ltrRoot();
	R=rho;
	delete &rho;
	computeGaussianCoefficients();
}



void
GPR::
computeGaussianCoefficients()
{
   for(int k=0;k<=N;k++) a[k]=EA(k,y);
}



//-----------THE COEFFICIENTS-----------------------------


// empirical coefficients
Real
GPR::
EC(int k, RealArray1D& w)       // w = y,yp is the data array
{
   assert((0<=k)&&(k<=N));

	Real sum=0.0;
	for(int j=0;j<=n;j++) sum+=w[j]*psi(k,j);
	return sum/(n+1);
}



// Gaussian coefficients
Real
GPR::
EA(int k, RealArray1D& w)         // w = y,yp is the data array
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

  // compute the Z_j, 0<=j<=n, gpr-notes (25), p. 11 with the
  // (unconditional) means E(E_j) added on the right hand side.
	RealArray1D Z(n+1);
	Z[0]=(w[0]-mu_E[0])/R(0,0);
	for(int j=1;j<=n;j++){

		Real sum=0.0;
		for(int k=0;k<j;k++) sum+=R(j,k)*Z[k];

		Z[j]=(w[j]-mu_E[j]-sum)/R(j,j);
	}

	// compute a_k=E[A_k|data], gpr-notes formula 25, p11 with the unconditional 
  // means added on the right, note A_k has (unconditional) mean mu[k]
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
   RealArray1D d=getCoefficients();
   Real f_q=0.0;
   for(int k=0;k<=q;k++) f_q+=d[k]*P[k];
   return f_q;
}



void
GPR::
expansionData(int q)
{
   assert(q<=N);       
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
      RealArray1D d=getCoefficients();
      Real sum=0.0;
      for(int i=0;i<=q;i++){ sum+=d[i]*psi[i]; fout << sum << "  "; }
      fout << endl;
   }
   fout.close();
   cout << endl << endl << endl
        << "Coefficients: " << getCoefficients() << endl;
}

			
//--------------------------TESTS------------------------------------		


Real
GPR::
min_R_jj()
{
   Real min=10000000.0;
   for(int j=0;j<=n;j++) if(R(j,j)<min) min=R(j,j);
   return min;
}


Real
GPR::
max_R_jj()
{
   Real max=0.0;
   for(int j=0;j<=n;j++) if(R(j,j)>max) max=R(j,j);
   return max;
}



Real
GPR::
conditionNumber_of_R()
{
   Real max=max_R_jj(),
        min=min_R_jj();
   if(max==0.0) return -1.0;
   return max/min;
}


void
GPR::
conditioning()
{
   bool done=false;
   while(!done){
     
      cout << "Condition of the Cholesky root R: " << conditionNumber_of_R() << endl
           << "R_jj min=" << min_R_jj() << ", R_jj max=" << max_R_jj() << endl
           << "Improve conditioning (0/1), improve=";
      int go; cin>>go;
      if(go==1){

        cout << "We'll add sigma to the diagonal of (K(s_i,s_j)). sigma=";
        Real sigma; cin>>sigma;
        diagonal_K_add(sigma);
      } else { done=true; }
   } // end while
}




void
GPR::
diagonal_K_add(Real delta)
{
   for(int j=0;j<=n;j++) K(j,j)+=delta;
   R=K.ltrRoot();
   if(have_gaussian) computeGaussianCoefficients();
}       

        
      

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
	Real f=1.0/sqrt(float(m+1)); B*=f;

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




//-----------------SET UP---------------------------------


RealArray1D
GPR::
dataPoints(int m, bool random)
{
   Real h=2.0/m;
   RealArray1D t(m+1);
   for(int i=0;i<=m;i++) t[i]=-1.0+i*h;

   if(random){

        vector<Real> u(m+1);
        u[0]=-1.0; u[m]=1.0;
        for(int i=1;i<m;i++) u[i]=-1.0+2.0*Random::U01();
        sort(u.begin(),u.end());
        for(int i=0;i<=m;i++) t[i]=u[i];
   }

   return t;
}


       
GPR&
GPR::
setUp()
{
   int N,n;
   Real sigma;
   RegressionType rType;
   BasisFunctions* bFcns;

   int ctr;   // prior mean centered? (0,1,2)
   Real r;    // 2: mu=(r,r,...,r)

   // DATA
   RealFunction  f;
   cout << "Sample function f:" << endl
        << "f0(t)=sin(2t\pi)...............[0]" << endl
        << "f1(t)=5t*exp(-9t^2/2)..........[1]" << endl
        << "f2(t)=(1+t)cos(8t\pi)..........[2]" << endl
        << "f3(t)=f1(t)+t^2sin(11t\pi).....[3]" << endl
        << "f4(t)=1.6|t|sin(10t\pi)........[4]" << endl
        << "f5(t)=t^2+t^3..................[5]" << endl << endl
        << "Enter f=[0,1,2,3,4,5]=";
   int fnum; cin>>fnum;
   // *f is reference to pointer to sample function
   f=functionExample(fnum);
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

   int random;
   cout << "Data points random (random=1) or evenly spaced (random=0), random=";
   cin>>random;

   // data abscissas s_j
   bool random_spacing=false;
   if(random==1) random_spacing=true;
   RealArray1D t=dataPoints(n,random_spacing);

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
   GPR* gpr=new GPR(N,t,f,bFcns,sigma,rType);

   // SET PRIOR MEAN
   if(rType==GPR::GAUSSIAN){

      // get center of Gaussian prior
      cout << "\nMean of Gaussian prior P:" << endl
           << "origin....................[0]" << endl
           << "empirical coefficients....[1]" << endl
           << "vector (r,r,...,r)........[2]" << endl
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
		lout << "Function data exact."
         << "Function: f" << fnum
         << endl;


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



//----------------CROSS VALIDATION-----------------------------


// still very inefficient since Cholesky root and kernel matrix fully recomputed.
// will not be fixed since this method is useless
int
GPR::
leaveOneOutCV()
{
    RealArray1D error(N+1);
    RealArray1D t(n);
    RealArray1D w(n);
    for(int k=0;k<=n;k++){       // leave out s_k

        for(int j=0;j<k;j++) { t[j]=s[j]; w[j]=y[j]; }
        for(int j=k;j<n;j++) { t[j]=s[j+1]; w[j]=y[j+1]; }

        // allocate the validating GPR
        GPR vGpr(N,t,w,basis,regrType);  
        const RealArray1D& a=vGpr.getCoefficients();
        const RealMatrix&  psi=vGpr.get_psi();
        // compute the prediction error at s_k for all the f_q
        Real f_q=0;
        for(int q=0;q<=N;q++){

           f_q+=a[q]*psi(q,k);     // f_q(s_j)
           error[q]+=(y[k]-f_q)*(y[k]-f_q);
        }
    } // error vector is computed

    int m=0; Real err=100000000;
    for(int q=0;q<=N;q++) if(error[q]<err){ err=error[q]; m=q; }

cout << endl << endl << "Leave one out error vector: " << error;
cout << endl << endl << "Optimal q: " << m;
    return m;
}
    

int
GPR::
polluteAndPredictCV()
{
    RealArray1D error(N+1);    // error[q]=error(f_q)

    // error is average over 3 trials
    for(int d=0;d<3;d++){

       // pollute data with Gaussian noise of standard deviation sigma
       for(int j=0;j<=n;j++) yp[j]=y[j]+sigma_*Random::sTN();
       RealArray1D d=computePollutedCoefficients();
       // compute the prediction error at s_k for all the f_q
       for(int k=0;k<=n;k++){     // error at s_k

          Real f_q=0;
          for(int q=0;q<=N;q++){     // error of f_q

              f_q+=d[q]*psi(q,k);     // f_q(s_j)
              error[q]+=(y[k]-f_q)*(y[k]-f_q);
          }       
       } // error vector at noise level sigma is computed
    } // end sigma

    int m=0; Real err=100000000;
    // enforce geometric decay of error as degree q of expansion increases
    // and discard f_N since this has inexplicably low errors
    for(int q=0;q<N;q++){

       // reduce roughness penalty in case of exact data.
       Real k_q = (sigma_<EPS) ? 1.0 : basis->roughnessPenalty(q);
       if(error[q]<k_q*err){ err=error[q]; m=q; }
    }
    return m;
}


const RealArray1D&
GPR::
computePollutedCoefficients()
{
   if(sigma_<EPS) return getCoefficients();
   // empirical polluted coefficients
   if(regrType==EMPIRICAL){

      for(int k=0;k<=N;k++) bp[k]=EC(k,yp);
      return bp;
   }
   // Gaussian polluted coefficients
   for(int k=0;k<=N;k++) ap[k]=EA(k,yp);
   return ap;
}





GPR_END_NAMESPACE(Gaussian)







   
