

#include "CubicSpline.hpp"
#include "TridiagonalMatrix.hpp"
#include <cmath>


namespace spline {
	

CubicSpline::
CubicSpline(int n_Data, Vector t, Vector w, double yp_Left, double yp_Right):
nData_(n_Data),
x_(t),
h_(n_Data-1),
y_(w),
ypLeft_(yp_Left),
ypRight_(yp_Right),
g_(n_Data),
sp_(n_Data-1)
{
		 int l=nData()-1;
		 Vector D(nData());
		 Vector diag(nData());

		 h(0) = x(1)-x(0);
		 D[0] = (y(1)-y(0))/h(0)-ypLeft();
		 diag[0] = 2*h(0);

		 for(int i=1;i<l;i++){

		 		 h(i) = x(i+1)-x(i);
		 		 D[i] = (y(i+1)-y(i))/h(i) - (y(i)-y(i-1))/h(i-1);
                 diag[i] = 2*(h(i)+h(i-1));
		 }
         D[l] = ypRight()-(y(l)-y(l-1))/h(l-1);
		 diag[l] = 2*h(l-1);

		 TridiagonalMatrix A(nData(),diag,h_);
		 A*=1/6.0;
		 g_=A.solve_AXeqY(D);

		 for(int i=0;i<l;i++) 
		 	sp(i) = (y(i+1)-y(i))/h(i) + h(i)*(g(i)+2*g(i+1))/6;
}


CubicSpline::
~CubicSpline(void)
{}


// Assumed that x\in[x(0),x(nData)]
// find i such that x\in(x(i),x(i+1)]
int
CubicSpline::
find(double t)
{
    // binary search
	int left=0, right=nData()-1;
	// maintain right>left and x\in(x(left),x(right)]
	while(right-left>1){

		 int m = left + (right-left)/2;
		 if(t>x(m)) left=m; else right=m;
	}
    return left;
}


double
CubicSpline::
value(double t)
{
    int l = nData()-1;
	if (t<x(0)) return y(0)+ypLeft()*(t-x(0));
	if (t>x(l)) return y(l)+ypRight()*(t-x(l));

	int i = find(t);
	double r = (g(i+1)-g(i))/h(i),
		   z = t-x(i+1);
	return y(i+1)+z*(sp(i)+(z/2)*(g(i+1)+r*z/3));
}


double
CubicSpline::
diff(double t)
{
    int l = nData()-1;
	if (t<x(0)) return ypLeft();
	if (t>x(l)) return ypRight();

	int i = find(t);
	double r = (g(i+1)-g(i))/h(i),
		   z = t-x(i+1);
	return sp(i)+z*(g(i+1)+r*z/2);
}


double
CubicSpline::
diff2(double t)
{
    int l = nData()-1;
	if (t<x(0)) return 0;
	if (t>x(l)) return 0;

	int i = find(t);
	double r = (g(i+1)-g(i))/h(i),
		   z = t-x(i+1);
	return g(i)+r*z;
}




} // end namespace spline
