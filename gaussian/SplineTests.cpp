


#include "SplineTests.hpp"
#include "TridiagonalMatrix.hpp"
#include "DoubleMatrix.hpp"
#include "CubicSpline.hpp"
#include "RNG.hpp"
#include <math.h>


namespace spline {


void
testTriangularSolver()
{
   cout << "Testing triangular solver..." << "\n\n";
   for(int k=0;k<20;k++){
  
      double error=1E-7;
      int dim = k+10;
      Vector Y(dim);
      Vector diag(dim);
      Vector subDiag(dim);  // extra element for convencience 

      for(int i=0;i<dim;i++){

		  Y[i]=nextUniform(-8,8);
          subDiag[i]= nextUniform(-2,2);
          diag[i]=nextUniform(6,12);         // diagonally dominant
      }
      LowerTriangularBidiagonalMatrix Q(dim,diag,subDiag);
      Vector X=Q.solve_QXeqY(Y);
      Vector Z=Q.rightMult(X);

      for(int i=0;i<dim;i++)
      if(fabs(Y[i]-Z[i]) > error) error = fabs(Y[i]-Z[i]);
      
      if (error>1E-7) cout << "Test failed, max error = " << error << endl;
      else cout << "Test passed." << endl;
   }
}



void
testTriangularSolver1()
{
   cout << "Testing transposed triangular solver..." << "\n\n";
   for(int k=0;k<20;k++){
  
      double error=1E-7;
      int dim = k+10;
      Vector Y(dim);
      Vector diag(dim);
      Vector subDiag(dim);  // extra element for convencience 

      for(int i=0;i<dim;i++){

		 		  Y[i]=nextUniform(-8,8);
         subDiag[i]= nextUniform(-2,2);
         diag[i]=nextUniform(6,12);         // diagonally dominant
      }
      LowerTriangularBidiagonalMatrix Q(dim,diag,subDiag);
      Vector X=Q.solve_QtXeqY(Y);
		   Vector Z=Q.qt_times_X(X);

      for(int i=0;i<dim;i++)
      if(fabs(Y[i]-Z[i]) > error) error = fabs(Y[i]-Z[i]);
      
      if (error>1E-7) cout << "Test failed, max error = " << error << endl;
      else cout << "Test passed." << endl;
   }
}



void
testSolver()
{
   srand(137);
   cout << "Testing bidiagonal solver..." << "\n\n";
   for(int k=0;k<20;k++){
  
      double error=1E-7;
      int dim = k+10;
      Vector Y(dim);
      Vector diag(dim);
      Vector subDiag(dim);  // extra element for convencience 

      for(int i=0;i<dim;i++){

		 Y[i]=nextUniform(-8,8);
         subDiag[i]= nextUniform(-2,2);
         diag[i]=nextUniform(6,12);         // diagonally dominant
      }
      TridiagonalMatrix A(dim,diag,subDiag);
      Vector X=A.solve_AXeqY(Y);
      Vector Z=A.rightMult(X);

      for(int i=0;i<dim;i++)
      if(fabs(Y[i]-Z[i]) > error) error = fabs(Y[i]-Z[i]);
      
      if (error>1E-7) cout << "Test failed, max error = " << error << endl;
      else cout << "Test passed." << endl;
   }
}



void
testGeneralTriangularSolver()
{
   cout << "Testing general solver..." << "\n\n";
   for(int k=0;k<20;k++){
  
      double error=1E-7;
      int dim = k+20;
      Vector Y(dim);
      LowerTriangularMatrix Q(dim);
	  double rho = 0.95;

      for(int i=0;i<dim;i++){

		 Y[i]=nextUniform(-8,8);
		 for(int j=0;j<=i;j++) Q(i,j)=exp(-rho*(i-j));   
      }
      Vector X=Q.solve_QXeqY(Y);
      Vector Z=Q.rightMult(X);

      for(int i=0;i<dim;i++)
      if(fabs(Y[i]-Z[i]) > error) error = fabs(Y[i]-Z[i]);
      
      if (error>1E-7) cout << "Test failed, max error = " << error << endl;
      else cout << "Test passed." << endl;
   }
}




void
testGeneralTriangularSolver1()
{
   cout << "Testing general transposed triangular solver..." 
	    << "\n\n";
   for(int k=0;k<20;k++){
  
      double error=1E-7;
      int dim = k+20;
      Vector Y(dim);
      LowerTriangularMatrix Q(dim);
	  double rho = 0.95;

      for(int i=0;i<dim;i++){

		 Y[i]=nextUniform(-8,8);
		 for(int j=0;j<=i;j++) Q(i,j)=exp(-rho*(i-j));   
      }
      Vector X=Q.solve_QtXeqY(Y);
      Vector Z=Q.transposeRightMult(X);

      for(int i=0;i<dim;i++)
      if(fabs(Y[i]-Z[i]) > error) error = fabs(Y[i]-Z[i]);
      
      if (error>1E-7) cout << "Test failed, max error = " << error << endl;
      else cout << "Test passed." << endl;
   }
}



void
testGeneralSolver()
{
   cout << "Testing general solver..." << "\n\n";
   for(int k=0;k<20;k++){
  
      double error=1E-7;
      int dim = k+20;
      Vector Y(dim);
      SymmetricMatrix A(dim);
	  double rho = 0.95;

      for(int i=0;i<dim;i++){

		 Y[i]=nextUniform(-8,8);
		 for(int j=0;j<=i;j++) A(i,j)=exp(-rho*(i-j));   
      }
      Vector X=A.solve_AXeqY(Y);
      Vector Z=A.rightMult(X);

      for(int i=0;i<dim;i++)
      if(fabs(Y[i]-Z[i]) > error) error = fabs(Y[i]-Z[i]);
      
      if (error>1E-7) cout << "Test failed, max error = " << error << endl;
      else cout << "Test passed." << endl;
   }
}






void
testCholeskyFactorisation()
{
   cout << "Testing tridiagonal Cholesky factorisation..." << "\n\n";
   for(int k=0;k<20;k++){

      int dim = k+10;
      Vector Y(dim);
      Vector diag(dim);
      Vector subDiag(dim);  // extra element for convencience 

      for(int i=0;i<dim;i++){

		 Y[i]=nextUniform(-8,8);
         subDiag[i]= nextUniform(-2,2);
         diag[i]=nextUniform(6,12);         // diagonally dominant
      }
      TridiagonalMatrix A(dim,diag,subDiag);
	  LowerTriangularBidiagonalMatrix Q=A.choleskyRoot();

	  double error = 1E-7;
	  for(int i=1;i<dim;i++){

		  double dg_i=Q.qqt().diag(i);
		  double sdg_i=Q.qqt().subDiag(i);
		  if(fabs(dg_i-A.diag(i)) > error) error = fabs(dg_i-A.diag(i));
		  if(fabs(sdg_i-A.subDiag(i)) > error) error = fabs(sdg_i-A.subDiag(i));

	  }
      if (error>1E-7) cout << "Test failed, max error = " << error << endl;
      else cout << "Test passed." << endl;
   }
}




void
testGeneralCholeskyFactorisation()
{
   cout << "Testing Cholesky factorisation..." << "\n\n";
   for(int k=0;k<20;k++){
  
      double error=1E-7;
      int dim = k+20;
      Vector Y(dim);
      SymmetricMatrix A(dim);
	  double rho = 0.95;

      for(int i=0;i<dim;i++){

		 Y[i]=nextUniform(-8,8);
		 for(int j=0;j<=i;j++) A(i,j)=exp(-rho*(i-j));   
      }
	  
 	  LowerTriangularMatrix Q = A.choleskyRoot();
	  SymmetricMatrix B = Q.qqt();
	 
	  for(int i=0;i<dim;i++)
	  for(int j=0;j<=i;j++)
		  if(fabs(A(i,j)-B(i,j))>error) error = fabs(A(i,j)-B(i,j));
	  if (error>1E-7) cout << "Test failed, max error = " << error << endl;
      else cout << "Test passed." << endl;
   }		  
}


double
nextUniform(double a, double b)
{
		 return a+(b-a)*rand()/RAND_MAX;
}

void 
testRand(double a, double b)
{
    for(int i=0;i<100;i++) cout << nextUniform(a,b) << endl;
}


void
testSpline()
{
		 int nData = 5;
		 
		 Vector X(nData);
		 Vector Y(nData);
		 
		 X[0]=-2; Y[0]=0.5;
         X[1]=-1.5; Y[1]=1;		 
		 X[2]=0.0; Y[2]=0.0;
         X[3]=1.2; Y[3]=2.5;		 
		 X[4]=2.0; Y[4]=0.5;
		 
		 double ypLeft=0;
		 double ypRight=0;

		 CubicSpline s(nData,X,Y,ypLeft,ypRight);
		 ofstream fout("spline.txt");
		 //ofstream fout_diff("spline_diff.txt");
		 //ofstream fout_diff2("spline_diff2.txt");

         // print: triple the number of data points.
		 // every third point is a knot.
		 double t=-2.0;
		 double dt=0.01;
		 
		 while(t<2.1){

			 fout << t 
			 	  << "; " << s(t) 
			      << "; " << s.diff(t)
			 	  << "; " << s.diff2(t)
			      << endl;
		 	 //fout_diff << t << "; " << s.diff(t) << endl;
		 	 //fout_diff2 << t << "; " << s.diff2(t) << endl;
		 	 t+=dt;
		 }
		 fout.close();
		 //fout_diff.close();
		 //fout_diff2.close();
}





} // end namespace spline
