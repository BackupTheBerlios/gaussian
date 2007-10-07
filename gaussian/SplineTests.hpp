

#include "TridiagonalMatrix.hpp"


namespace spline {

/** Test solver QX=Y, Q lower triangular and bidiagonal.*/
void testTriangularSolver();
/** Test solver Q'X=Y, Q lower triangular and bidiagonal.*/
void testTriangularSolver1();	
/** Test solver QX=Y, Q lower triangular. */
void testGeneralTriangularSolver();
/** Test solver Q'X=Y, Q lower triangular.*/
void testGeneralTriangularSolver1();
/** Test solver AX=Y, A tridiagonal.*/
void testSolver();
/** Test solver AX=Y, A symmetric.*/
void testGeneralSolver();
/** Cholesky factorization for symmetric tridiagonal matrices.*/
void testCholeskyFactorisation();
/** Cholesky factorization for symmetric matrices.*/
void testGeneralCholeskyFactorisation();

double nextUniform(double a, double b);
void testRand(double a, double b);

void testSpline();
void testQPInterpolator();


} // end namespace spline
