
#include "TridiagonalMatrix.hpp"
#include <cmath>


namespace spline {



LowerTriangularBidiagonalMatrix::
LowerTriangularBidiagonalMatrix(int dim)
:
dim_(dim)
{
    entries_.resize(dim);
    for(int i=0;i<dim;i++) entries_[i].resize(2);
}


LowerTriangularBidiagonalMatrix::
LowerTriangularBidiagonalMatrix(int dim, Vector diag, Vector subDiag)
:
dim_(dim)
{
    entries_.resize(dim);
    // no subdiagonal element
    entries_[0].resize(2);
    entries_[0][0]=0.0;
    entries_[0][1]=diag[0];
    // subdiagonal has length dim-1
    for(int i=1;i<dim;i++){

         entries_[i].resize(2);
         entries_[i][0]=subDiag[i-1];
         entries_[i][1]=diag[i];
    }
}




double 
LowerTriangularBidiagonalMatrix::
operator()(int i,int j) const
{
    if(i==j) return entries_[i][1];
    if((i==j-1)||(i==j+1)) return entries_[i][0];

    return 0;
}


Vector 
LowerTriangularBidiagonalMatrix::
solve_QXeqY(Vector Y) const
{
    Vector X(dim());
    X[0]=Y[0]/diag(0);
    for(int i=1;i<dim();i++) X[i]=(Y[i]-subDiag(i)*X[i-1])/diag(i);

    return X;
}


Vector 
LowerTriangularBidiagonalMatrix::
solve_QtXeqY(Vector Y) const
{
    Vector X(dim());
	int l = dim()-1; X[l]=Y[l]/diag(l);
    for(int i=l-1;i>=0;i--) X[i]=(Y[i]-subDiag(i+1)*X[i+1])/diag(i);

    return X;
}



Vector 
LowerTriangularBidiagonalMatrix::
rightMult(Vector X) const
{
	Vector Y(dim());
	Y[0]=diag(0)*X[0];
	for(int i=1;i<dim();i++) 
	Y[i]=subDiag(i)*X[i-1]+diag(i)*X[i];

    return Y;
}



Vector 
LowerTriangularBidiagonalMatrix::
qt_times_X(Vector X)
{
	Vector Y(dim());
	for(int i=0;i<dim()-1;i++) 
	Y[i]=diag(i)*X[i]+subDiag(i+1)*X[i+1];

	int k=dim()-1;
    Y[k]=diag(k)*X[k];

	return Y;
}



TridiagonalMatrix 
LowerTriangularBidiagonalMatrix::
qqt()
{
	TridiagonalMatrix A(dim());
	A.diag(0) = diag(0)*diag(0);
	A.subDiag(0)=diag(0)*subDiag(1);

	for(int i=1;i<dim();i++){
 
		A.diag(i) = diag(i)*diag(i)+subDiag(i)*subDiag(i);
		A.subDiag(i)=diag(i-1)*subDiag(i);
	}
	return A;
}





//--------------------------- tridiagonalmatrix ----------


TridiagonalMatrix::
TridiagonalMatrix(int dim)
:
dim_(dim)
{
    entries_.resize(dim);
    for(int i=0;i<dim;i++) entries_[i].resize(2);
}



TridiagonalMatrix::
TridiagonalMatrix(int dim, Vector diag, Vector subDiag)
:
dim_(dim)
{
    entries_.resize(dim);
    // no subdiagonal element
    entries_[0].resize(2);
    entries_[0][0]=0.0;
    entries_[0][1]=diag[0];
    // subdiagonal has length dim-1
    for(int i=1;i<dim;i++){

         entries_[i].resize(2);
         entries_[i][0]=subDiag[i-1];
         entries_[i][1]=diag[i];
    }
}



double 
TridiagonalMatrix::
operator()(int i,int j) const
{
    if(i==j) return entries_[i][1];
    if((i==j-1)||(i==j+1)) return entries_[i][0];

    return 0;
}



LowerTriangularBidiagonalMatrix 
TridiagonalMatrix::
choleskyRoot() const
{
	LowerTriangularBidiagonalMatrix Q(dim());
	Q.diag(0) = sqrt(diag(0));
	Q.subDiag(0) = 0.0;
	for(int i=1;i<dim();i++){
		 		 
		Q.subDiag(i)=subDiag(i)/Q.diag(i-1);
		Q.diag(i)=sqrt(diag(i)-Q.subDiag(i)*Q.subDiag(i));
	}
	return Q;
}



Vector 
TridiagonalMatrix::
solve_AXeqY(Vector Y) const
{
   LowerTriangularBidiagonalMatrix Q=choleskyRoot();
   Vector Z=Q.solve_QXeqY(Y);
   return Q.solve_QtXeqY(Z);
}



Vector 
TridiagonalMatrix::
rightMult(Vector X) const
{
	Vector Y(dim());
	Y[0]=diag(0)*X[0]+subDiag(1)*X[1];
	for(int i=1;i<dim()-1;i++) 
		Y[i]=subDiag(i)*X[i-1]+diag(i)*X[i]+subDiag(i+1)*X[i+1];
	int k = dim()-1;
	Y[k]=subDiag(k)*X[k-1]+diag(k)*X[k];

	return Y;
}




bool 
TridiagonalMatrix::
testAXeqY(Vector X, Vector Y)
{
	double error=1E-7;
	Vector Z=rightMult(X);

	for(int i=0;i<dim();i++)
	if(fabs(Y[i]-Z[i]) > error) error = fabs(Y[i]-Z[i]);
      
	if (error>1E-7) return false;
	else return true;
}



void 
TridiagonalMatrix::
operator*=(double t)
{
	for(int i=0;i<dim();i++){

		entries_[i][0]*=t;
		entries_[i][1]*=t;
	}
}




} // end namespace spline
