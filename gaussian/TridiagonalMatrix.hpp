

#ifndef SPLINETRIDIAGONALMATRIX_H
#define SPLINETRIDIAGONALMATRIX_H

#include "Typedefs.hpp"



namespace spline {


class TridiagonalMatrix;  // forward declaration
/**
 * @author mjhmeyer
 * Lower triangular matrix for which only the diagonal and subdiagonal are nonzero.
 * Only these two diagonals are stored. 
 */
class LowerTriangularBidiagonalMatrix {

public:

explicit LowerTriangularBidiagonalMatrix(int dim);
LowerTriangularBidiagonalMatrix(int dim, Vector diag, Vector subDiag);
//~TridiagonalMatrix(){};

/** Entries. Slow. Use {@link diag}, {@link subDiag} instead.
 */
double operator()(int i,int j) const;
/** QX, where Q=this.
 */
Vector rightMult(Vector X) const;
/** Q'X, where Q=this.
 */
Vector qt_times_X(Vector X);
int dim() const { return dim_; }
/** diagonal element this[i][i], read-write access.
 */
double& diag(int i){ return entries_[i][1]; }
/** subdiagonal element this[i-1][i] (=0, if i=0), read-write access.
 */
double& subDiag(int i){ return entries_[i][0]; }
/** diagonal element this[i][i], read-access.
 */
double diag(int i) const { return entries_[i][1]; }
/** subdiagonal element this[i-1][i] (=0, if i=0), read-access.
 */
double subDiag(int i) const { return entries_[i][0]; }
/** Solves the equation \f$QX=Y\f$, where Q=this.
 *  Does not check for nonsingularity.
 * @returns the solution X.
 */
Vector solve_QXeqY(Vector Y) const;
/** Solves the equation \f$Q'X=Y\f$, where Q=this.
 *  Does not check for nonsingularity of Q.
 * @returns the solution X.
 */
Vector solve_QtXeqY(Vector Y) const;

TridiagonalMatrix qqt();


private:

int dim_;
// 2 by dim matrix
Matrix entries_;    // [j][0]=matrix(j-1,j), [j][1]=diag(j)

};




/**
 * @author mjhmeyer
 * Symmetric square tridiagonal matrix.
 * Only the diagonal and the subdiagonal are stored.
 */
class TridiagonalMatrix {

public:

explicit TridiagonalMatrix(int dim);
TridiagonalMatrix(int dim, Vector diag, Vector subDiag);
~TridiagonalMatrix(){};

/** Entries. Slow. Use {@link diag}, {@link subDiag} instead.
 */
double operator()(int i,int j) const;
/** Scalar multiplication of this with t. */
void operator*=(double t);
/** Multiply this with X from the right */
Vector rightMult(Vector X) const;
int dim() const { return dim_; }
/** diagonal element this[i][i], read-write access.
 */
double& diag(int i){ return entries_[i][1]; }
/** subdiagonal element this[i-1][i] (=0, if i=0), read-write access.
 */
double& subDiag(int i){ return entries_[i][0]; }
/** diagonal element this[i][i], read only access.
 */
double diag(int i) const { return entries_[i][1]; }
/** subdiagonal element this[i-1][i] (=0, if i=0), read only access.
 */
double subDiag(int i) const { return entries_[i][0]; }
/** Returns a lower triangular matix Q such that QQ'=this. 
 *  Only the diagonal and first subdiagonal are nonzero.
 */
LowerTriangularBidiagonalMatrix choleskyRoot() const;


/** Solves the equation \f$AX=Y\f$, where A=this using the
 *  Cholesky factorization of A.
 *  Does not check for nonsingularity of A.
 * @returns the solution X.
 */
Vector solve_AXeqY(Vector Y) const;

/** tests if AX=Y, where A=this;
 */
bool testAXeqY(Vector X, Vector Y);


private:

int dim_;
// 2 by dim matrix
Matrix entries_;    // [j][0]=matrix(j-1,j), [j][1]=diag(j)

};




} // end namespace spline




#endif
