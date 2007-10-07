/***************************************************************************
 *            Matrix.cpp
 *
 *  Sun Oct  7 20:12:57 2007
 *  Copyright  2007  mjhmeyer@yahoo.com
 *  
 ****************************************************************************/
 

#include "DoubleMatrix.h"
#include <cmath>


SymmetricMatrix::
SymmetricMatrix(int n)
:
dim_(n)
{
	A.resize(n);
	for(int i=0;i<n;i++) A[i].resize(i+1);	
}


double& 
SymmetricMatrix::
operator()(int i, int j)
{
	if(i<=j) return A[i][j];
	return A[j][i];	
}


LowerTriangularMatrix 
SymmetricMatrix::
choleskyRoot()
{
	LowerTriangularMatrix Q(dim());
	double sum = 0.0;
	
    Q(0,0)=sqrt(A[0][0]);
    for(int i=0;i<dim();i++){

		// Qij: row(j).row(i) = Aij
       	for(int j=0;j<i;j++){

			sum = 0.0;
		   	for(int k=0;k<j;k++) sum+=Q(j,k)*Q(i,k);
			Q(i,j)=(A[i][j]-sum)/Q(j,j);		   
	   	}
	 	// Qii: row(i).row(i) = Aii
		sum = 0.0;
		for(int k=0;k<i;k++) sum+=Q(i,k)*Q(i,k);
		Q(i,i) = sqrt(A[i][i]-sum);
	}		
	return Q;	
}


Vector 
SymmetricMatrix::
rightMult(Vector X)
{
	Vector Y(dim());
    for(int i=0;i<dim();i++){

		double sum = 0.0;
		for(int j=0;j<dim();j++) sum += A[i][j]*X[j];
		Y[i] = sum;
	}		
	return Y;	
}


Vector 
SymmetricMatrix::
solve_AXeqY(Vector Y)
{
	LowerTriangularMatrix Q = choleskyRoot();	
	Vector W = Q.solve_QXeqY(Y);
	return Q.solve_QtXeqY(W);
}



Vector 
LowerTriangularMatrix::
rightMult(Vector X)
{
	Vector Y(dim());
    for(int i=0;i<dim();i++){

		double sum = 0.0;
		for(int j=0;j<=i;j++) sum += Q[i][j]*X[j];
		Y[i] = sum;
	}		
	return Y;	
}


SymmetricMatrix 
LowerTriangularMatrix::
qqt()
{
	SymmetricMatrix A(dim());
	for(int i=0;i<dim();i++)
	for(int j=0;j<=i;j++){
		
		double sum = 0.0;
		for(int k=0;k<=j;k++) sum += Q[i][k]*Q[j][k];
		A(i,j)=sum;
	}
	return A;
}


Vector 
LowerTriangularMatrix::
solve_QXeqY(Vector Y)
{
	Vector X(dim());
	for(int i=0;i<dim();i++){
		
		double sum=0.0;
		for(int j=0;j<i;j++) sum += Q[i][j]*X[j];
		X[i] = (Y[i]-sum)/Q[i][i];			
	}
	return X;
}



Vector 
LowerTriangularMatrix::
solve_QtXeqY(Vector Y)
{
	Vector X(dim());
	for(int i=dim()-1;i>=0;i--){
		
		double sum=0.0;
		for(int j=i+1;j<dim();j++) sum += Q[i][j]*X[j];
		X[i] = (Y[i]-sum)/Q[i][i];			
	}
	return X;
}
