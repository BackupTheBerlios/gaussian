 /***************************************************************************
 *            gpr.h
 *
 *  Sat May 29 16:37:19 2004
 *  Copyright  2004  cpp
 *  cpp@linux
 ****************************************************************************/


#ifndef gpr_gpr_h
#define gpr_gpr_h


#include "Matrix.h"
#include "BasisFunctions.h"
#include "FunctionExamples.h"
#include <string>
using std::string;


GPR_BEGIN_NAMESPACE(Gaussian)


class Functional;
class LinearFunctional;

/** Empirical and Gaussian process regression on [-1,+1] using expansion in
 *  basis functions. See gprs-1.2.ps ("gpr-notes").
 **/
class GPR {
	
public:

   /** Flag, center of Gaussian prior (origin or empirical coefficients). */
   enum RegressionType {EMPIRICAL, GAUSSIAN};

   /** Maximal expansion length: \f$\psi_0,\psi_1,\dots,\psi_N\f$. */
   int get_N() const { return N; }

   /** Parameter n = number of data points - 1. */
   int get_n() const { return n; }

   /** The points \f$s_j\f$. */ 
   const RealArray1D& get_s() const { return s; }

   /** The values \f$y_j=f(s_j)+noise\f$. */
   const RealArray1D& getFunctionData() const { return y; }

   /** Current regression type (EMPIRICAL, GAUSSIAN). */
   RegressionType getRegressionType() const { return regrType; }
   
   /** Set the type of regression (EMPIRICAL, GAUSSIAN).
    *  Don't set it to GAUSSIAN on very large data sets (>1000 points).
    *  The matrix omputations are too slow in this case.
    */
   void setRegressionType(RegressionType rt);

   /** The matrix \f$(\psi_i(s_j))\f$. */
   const RealMatrix& get_psi() const { return psi; }

   /** The vector of empirical coefficients. */
   const RealArray1D& getEmpiricalCoefficients() const { return b; }

   /** The vector of Gaussian coefficients. */
   const RealArray1D& getGaussianCoefficients() const { return a; }

   /** The vector of coefficients used by the current regression. */
   const RealArray1D& getCoefficients() const
   { if(regrType==EMPIRICAL) return b; return a; }

   /** Is the Gaussian machinery enabled?
    */
   bool haveGaussian() const { return have_gaussian; }

   /** The mean of the Gaussian prior P encoded as the vector
    *  \f$\mu_k=E[A_k], k\leq N,\f$ of the unconditional 
    *  expectations of the coefficient functionals.
    */
   const RealArray1D& getPriorMean() const { return mu; }

   /** Set the mean \f$\mu_k=E[A_k]\f$ of the Gaussian prior P equal to 
    *  mean. Does nothing if regression type is empirical.
    */
   void setPriorMean(const RealArray1D& mean);

   /** Set the mean of the Gaussian prior P equal to mean.
    *  Does nothing if regression type is empirical.
    */
   void setPriorMeanToOrigin();

   /** Set the mean of the Gaussian prior P equal to mean.
    *  Does nothing if regression type is empirical.
    */
   void setPriorMeanToEmpiricalCoefficients();

   /** Set the basis functions (must recompute all data structures).*/
   const BasisFunctions* getBasisFunctions() const { return basis; };

   /** Set the basis functions (must recompute all data structures).*/
   void setBasisFunctions(BasisFunctions* bFcns);

   /** Set the function data array.*/
   void setFunctionData(RealArray1D& w);

   /** Set the function underlying the data. This fills the function data
    *  array with <i>exact</i> data. Add noise if desired.
    */
   void setFunction(RealFunction g);

   /** Adds independent gaussian noise of standard deviation sigma
    *  to the function data (recomputes all coefficients).
    */
   void addNoise(Real sigma);


//---------------CONSTRUCTORS---------------------------

   /** The regression type is set to GAUSSIAN and the
    *  prior mean to the origin. Use public methods to change this
    *  setup. If the data set is very large (<1000 points) we may not
    *  want to enable the Gaussian machinery (large matrices) and try to get
    *  by with empirical regression.
    *
    *  @param dmax maximum degree N of expansion.
    *  @param t data abscissas \f$t_j,\ 0<=j<=n\f$.
    *  @param w function data \f$w_j=f(t_j),\ 0<=j<=n\f$.
    *  @param bFcns object containing the basis functions.
    *  @param sigma standard deviation of noise in the data.
    *  @param rt type of regression (EMPIRICAL, GAUSSIAN),
    *  EMPIRICAL does not initialize the Gaussian data structures.
    **/
   GPR(int dmax, RealArray1D& t, RealArray1D& w, BasisFunctions* bFcns,
       Real sigma=0.2, RegressionType rt=EMPIRICAL);

   /** As {@link GPR(int,RealArray1D&,RealArray1D&,BasisFunctions*,RegressionType)}
    *  except that the function generating the data is known. This is for
    *  tests of the algorithm.
    *
    *  @param g function generating the data.
    *  @param sigma standard deviation of data noise.
    **/
   GPR(int dmax, RealArray1D& t, RealFunction g, BasisFunctions* bFcns, 
       Real sigma=0.2, RegressionType rt=EMPIRICAL);


   /** Sets up empirical regression with clean data from g.
    *
    * @param random data points \f$s_j\f$ evenly spaced or random.
    * @param M M+1 basis functions.
    * @param m m+1 data points.
    */
   GPR(int M, int m, RealFunction g, BasisFunctions* bFcns, bool random);


//--------------------SETUP-------------------------------------------------


   /** Sets up a Gaussian process regressor {@link GPR} with function examples
    *  via a dialogue with the user. User must delete the GPR object.
    */
   static GPR& setUp();

   /** Array of points \f$t_j\in[-1,+1],\ 0\leq j\leq m\f$ with
    *  \f$t_0=-1,\ t_m=+1\f$.
    *  Points are in increasing order and equidistant or uniformly random.
    *  Use: setup of a GPR object.
    */
   static RealArray1D dataPoints(int m, bool random=false);


//--------------------EXPANSIONS--------------------------------------------


   /** The sequence \f$(\psi_0(t),\psi_1(t),\dots,\psi_m(t))\f$ of the
    *  basis functions evaluated at t.
    */
   RealArray1D basisFunctionValues(Real t, int m)
   { return basis->values(t,m); }


	/** The expansion \f$f_q(t)\f$ using Gaussian or empirical coefficients \f$a_k\f$
	 *  depending on the current state of <code>this</code>.
	 **/
	Real expansion(Real t, int q);

  /** Writes data files for the expansions expansions \f$f_0,f_1,\dots,f_q\f$
   *   of f in Legendre polynomials on [-1,+1] with coefficients computed by
   *  Gaussian regression for potting with gnuplot.
   *
   *  The expansions are evaluated at 801 evenly spaced points \f$t_j\in[-1+1]\f$
   *  and written to a file "ExpansionData.txt" in gnuplot data format.
   *  The first column contains the points \f$t_j\f$, the next column the true
   *  function values \f$f(t_j)\f$ (if the function f generating the data is known
   *  otherwise skipped) and the subsequent columns the expansions
   *  \f$f_0(t_j),f_1(t_j),\dots,f_q(t_j)\f$.
   *
   *  See "doc/gnuplot_Readme.html" for instructions how to plot such data
   *  with gnuplot.
   *
   *  The data points are written to the file "FunctionData.txt". The first column
   *  contains the points \f$s_j\f$ and the second column the data points \f$y_j\f$.
   *
   **/
   void expansionData(int q);


//-----------------CROSS VALIDATION---------------------------------


   /** Leave one out cross validation.
    *  @returns index q of \$f_q\f$ with minimal error.
    */
   int leaveOneOutCV();

   /** The data are polluted with independent Gaussian noise (\f$\sigma=0.2\f$),
    *  the regressors computed and we then test which regressor best predicts the 
    *  original data, see gpr-notes, section 5.2, p. 32.
    *
    *  @returns index q of \$f_q\f$ with minimal error.
    */
   int polluteAndPredictCV();



//------------------FUNCTIONAL ESTIMATION----------------------------


    /** Estimate the value \f$E[\thinspace L\thinspace|\thinspace data\thinspace]\f$ 
     *  as in gpr-notes, section 6 (pp. 33-35). L has to be based on the
     *  GPR <code>this</code> (L-constructor).
     *
     *  @returns \f$E[\thinspace L\thinspace|\thinspace data\thinspace]\f$.
     */
    Real estimateFunctional(Functional& L);


    /** Estimate the value \f$E[\thinspace L\thinspace|\thinspace data\thinspace]\f$ 
     *  as in gpr-notes, section 6, equation (40), p. 34. This assumes that the
     *  coefficients \f$a_k=E[\thinspace A_k\thinspace|\thinspace data\thinspace]\f$ 
     *  have already been computed. L has to be based on the GPR <code>this</code>.
     *
     *  Current implementation uses the regressor \f$f_N\f$. This is only
     *  useful with exact data.
     *
     *  @returns \f$E[\thinspace L\thinspace|\thinspace data\thinspace]\f$.
     */
    Real estimateLinearFunctional(LinearFunctional& L);


    /** Monte Carlo integral over [-1,+1] of the underlying function f computed 
     *  from the function data. Use this to compare performance with Bayesian
     *  Monte Carlo.
     */
    Real monteCarloIntegral();

   


//------------------TESTS--------------------------------------------

   /** Reports the condition number of the Cholesky root R of the kernel matrix
    *  \f$D=(K(s_i,s_j))\f$ and asks the user if the conditioning should be improved
    *  by adding $\sigma>0$ to the diagonal elements of $D$.
    */
   void conditioning();

   /** Condition number of the Cholesky root R of the kernel matrix
    *  \f$(K(s_i,s_j))\f$.
    */
   Real conditionNumber_of_R();
   
   /** Allocates a directory <code>bases</code> and prints a plot of the
    *  first q basis functions for each of the bases which are implemented.
    **/
   static void printBasisFunctions(int q);


   /** Tests the basis functions \f$\psi_0,\dots,\psi_q\f$ for orthonormality
    *  using Monte Carlo integration on m random points in [-1,+1].
    *  Prints the matrix of inner products \f$(\psi_i,\psi_j)\f$.
    **/
   void orthoTest(int q, int m);

  

protected: 


   /** The lower triangular Cholesky root of the kernel matrix 
    *  \f$(K(s_i,s_j))\f$.
    */
   LTRRealMatrix& getR(){ return R; }
	
	

private:

   static Real EPS;         // standard deviation smaller than this considered zero

   bool have_gaussian;      // flag, Gaussian machinery available
   RegressionType regrType; // flag, type of regression (EMPIRICAL, GAUSSIAN).
   int N;                   // expansions cut off at psi_N
   int n;                   // points s_j, j=0,...,n.
   Real sigma_;             // noise level in data
   RealArray1D s;           // the points s_j, j<=n
   RealArray1D y;           // the values y_j, j<=n
   RealArray1D yp;          // polluted data y_j+z_j, j<=n for reverse prediction
   RealFunction f;          // function generating the data (NULL if unknown).
   BasisFunctions* basis;   // the basis functions psi_k
   string basis_name;       // name of basis
   RealMatrix psi;          // basis functions psi(k,j)=psi_k(s_j), 0<=k<=N, 0<=j<=n.
   UTRRealMatrix K;         // kernel matrix K(s_i,s_j), 0<=i,j<=n.
   LTRRealMatrix R;         // lower triangular Cholesky root of K.
   RealArray1D mu;          // center of the Gaussian prior: mu[k]=E[A_k].
   RealArray1D a;           // Gaussian coefficients
   RealArray1D b;           // empirical coefficients
   RealArray1D ap;          // Gaussian coefficients based on the polluted data yp
   RealArray1D bp;          // empirical coefficients based on the polluted data yp

   // Gaussian coefficient a_k=E(A_k) of psi_k based on data w (=y,yp),
   // we need the data array yp for "pollute and predict" cross validation.
   // see gpr-notes, section 5.2, p. 32.
   Real EA(int k, RealArray1D& w);
   // empirical coefficient of psi_k  based on data w (=y,yp):      
   Real EC(int k, RealArray1D& w);
         
   Real min_R_jj();
   Real max_R_jj();
   void diagonal_K_add(Real delta);
   
   void initBasis();                    // matrix of basis function values
   void computeEmpiricalCoefficients(); // write into coefficient array
   void initGaussian();                 // all data for Gaussian regression
   void computeGaussianCoefficients();  // write into coefficient array
   // the coefficients based on the polluted data yp
   const RealArray1D& computePollutedCoefficients();

}; // end GPR



GPR_END_NAMESPACE(Gaussian)

 

#endif
