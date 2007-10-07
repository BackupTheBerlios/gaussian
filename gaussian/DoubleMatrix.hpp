/***************************************************************************
 *            Matrix.h
 *
 *  Sun Oct  7 19:58:15 2007
 *  Copyright  2007  mjhmeyer@yahoo.com
 *  
 ****************************************************************************/

#ifndef ____Matrix_H_1949532____
#define ____Matrix_H_1949532____


/*! @file DoubleMatrix.hpp
 * Symmetric matrices of type double, Cholesky factorization and 
 * solution of systems of equations.
 */
 
#include "Typedefs.hpp"
 
class LowerTriangularMatrix;
 
 
class SymmetricMatrix {
 
public:
 
SymmetricMatrix(int n);
~SymmetricMatrix(){};
 
int dim(){ return dim_; }
LowerTriangularMatrix choleskyRoot();
Vector rightMult(Vector X);
Vector solve_AXeqY(Vector Y);
/** Read-write access, no bounds checking.*/    
double& operator()(int i, int j);
 
private:

int dim_;
Matrix A;

};




class LowerTriangularMatrix {
 
public:
 
LowerTriangularMatrix(int n);
~LowerTriangularMatrix(){}
 
int dim(){ return dim_; }
/* QX where Q=this. */
Vector rightMult(Vector X);
/* Q'X where Q=this. */
Vector transposeRightMult(Vector X);
/* QQ' where Q=this. */
SymmetricMatrix qqt();
/* Solves QX=Y, where Q=this, and returns the solution. */
Vector solve_QXeqY(Vector Y);
/* Solves Q'X=Y, where Q=this, and returns the solution. */
Vector solve_QtXeqY(Vector Y);
/** Read-write access.
 *  i,j must satisy i<=j, no bounds checking. 
 */    
double& operator()(int i, int j);
 
private:

int dim_; 
Matrix Q;
  
};


#endif
