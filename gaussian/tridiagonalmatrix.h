/***************************************************************************
 *   Copyright (C) 2007 by mjhmeyer   *
 *   mjhmeyer@yahoo.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef SPLINETRIDIAGONALMATRIX_H
#define SPLINETRIDIAGONALMATRIX_H

#include "Typedefs.h"

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
Vector operator*=(Vector X) const;
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
Vector operator*=(Vector X) const;
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


private:

int dim_;
// 2 by dim matrix
Matrix entries_;    // [j][0]=matrix(j-1,j), [j][1]=diag(j)

};

} // end namespace spline

#endif
