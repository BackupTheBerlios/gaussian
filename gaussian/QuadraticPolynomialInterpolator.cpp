/***************************************************************************
 *            QuadraticPolynomialInterpolator.cpp
 *
 *  Sun Oct  7 23:22:54 2007
 *  Copyright  2007  mjhmeyer@yahoo.com
 *  
 ****************************************************************************/


#include "QuadraticPolynomialInterpolator.hpp"
#include "DoubleMatrix.hpp"



QuadraticPolynomialInterpolator::
QuadraticPolynomialInterpolator(int nData, Vector x, Vector y)
{
	double 	sum1=0.0, 
			sum2=0.0, 
			sum3=0.0, 
			sum4=0.0;
	
	SymmetricMatrix A(3);
	for(int i=0;i<nData;i++){
		
		sum1+=x[i];
		sum2+=x[i]*x[i];
		sum3+=x[i]*x[i]*x[i];
		sum4+=x[i]*x[i]*x[i]*x[i];
	}
	A(0,0)=nData;
	A(0,1)=sum1;
	A(0,2)=sum2;
	A(1,1)=sum2;
	A(1,2)=sum3;
	A(2,2)=sum4;
	
	sum1=sum2=sum3=0.0;
	for(int i=0;i<nData;i++){
		
		sum1+=y[i];
		sum2+=x[i]*y[i];
		sum3+=x[i]*x[i]*y[i];
	}
	Vector W(3);
	W[0]=sum1; W[1]=sum2; W[2]=sum3;
	
	a.resize(3);
	a=A.solve_AXeqY(W);	
}




double 
QuadraticPolynomialInterpolator::
operator()(double t)
{
	return a[0]+t*(a[1]+t*a[2]);		
}



double 
QuadraticPolynomialInterpolator::
diff(double t)
{
	return a[1]+2*t*a[2];		
}



double 
QuadraticPolynomialInterpolator::
diff2(double t)
{
	return 2*a[2];		
}
