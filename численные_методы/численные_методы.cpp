using namespace std;
#include <iostream>
#include "функции.h"
#define Omega 1.
double f4(double x)
{
	return x * x * x * sqrt(4 + x * x);
}
double f31(double* x)
{
	return exp(x[0] * x[1]) + x[0] - 4;
}
double f32(double* x)
{
	return x[0]*x[0] -4*x[1] -1;
}
double f31_(double* x)
{
	return -exp(x[0] * x[1])  + 4;
}
double f32_(double* x)
{
	return (x[0] * x[0] - 1) / 4;
}
double df311(double* x)
{
	return x[1]*exp(x[0] * x[1]) + 1;
}
double df322(double* x)
{
	return 2 * x[0] ;
}
double df313(double* x)
{
	return x[1]*exp(x[0] * x[1]);
}
double df324(double* x)
{
	return - 4 ;
}
double f2(double x)
{
	return log10(2 * x + 1) - x * x * x + 1;
}
double f2_(double x)
{
	return sqrt(sqrt(sqrt((log10(2 * x + 1)  + 1)* (log10(2 * x + 1) + 1))));
}
double df2(double x)
{
	return 2 / (2 * x + 1) / log(10) - 3 * x * x;
}
double f1(double x)
{
	return x*cos(x);
}
int main()
{
    setlocale(LC_ALL, "Russian");
	/*int n = 4;
	double pi = 3.14;
	double* x = new double[n];
	double* y = new double[n];
	double X = pi / 4;
	x[0] = 0;
	x[1] = pi / 6;
	x[2] = pi / 3;
	x[3] = pi / 2;
	for (int i = 0; i < n; i++)
	{
		y[i] = f1(x[i]);
	}
	cout << inter_laplas(x, y, n, X) << endl;
	cout << inter_Newton(x, y, n, X) << endl;
	x[0] = 0;
	x[1] = pi / 6;
	x[2] = 5*pi / 12;
	x[3] = pi / 2;
	for (int i = 0; i < n; i++)
	{
		y[i] = f1(x[i]);
	}
	cout << inter_laplas(x, y, n, X) << endl;
	cout << inter_Newton(x, y, n, X) << endl;*/
	/*int n = 5;
	double* x = new double[n];
	double* y = new double[n];
	double X = 1.5;
	x[0] = 0;
	x[1] = 1;
	x[2] = 2;
	x[3] = 3;
	x[4] = 4;
	y[0] = 0;
	y[1] = 0.45345;
	y[2] = 0.5236;
	y[3] = 0;
	y[4] = -2.2672;
	cout << cubic_spline(x, y, n, X) << endl; */
	//int n = 5;
	//int m = 3;
	//double* x = new double[n];
	//double* y = new double[n];
	//double* res = new double[n];
	//double X = 1.5;
	//x[0] = -1;
	//x[1] = 0;
	//x[2] = 1;
	//x[3] = 2;
	//x[4] = 3;
	//y[0] = -0.86603;
	//y[1] = 0;
	//y[2] = 0.86603;
	//y[3] = 1;
	//y[4] = -4.3301;
	//MNK(x,y,res,n,m);
	//for (int i = 0; i < m; i++)
	//{
	//	cout << res[i] << "*x^" << i << endl;
	//}
	/*cout << m_hord(f2, -1, 0) << endl;
	cout << m_hord(f2, 1, 2) << endl;
	cout << m_newtona(f2, df2, 0) << endl;
	cout << m_pros_it(f2_, 1, 2);*/
	/*int n = 2;
	f* f1 = new f[n];
	f** df1 = new f*[n];
	double* res = new double[n];
	for (int i = 0; i < n; i++)
	{
		df1[i] = new f[n];
	}
	f1[0] = f31;
	f1[1] = f32;
	df1[0][0] = df311;
	df1[0][1] = df313;
	df1[1][0] = df322;
	df1[1][1] = df324;
	newton_sis(f1, df1, n, res);
	for (int i = 0; i < n; i++)
	{
		cout << res[i] << endl;
	}*/
	/*int n = 2;
	double* res = new double[n];
	f* f1 = new f[n];
	f1[0] = f31_;
	f1[1] = f32_;
	pr_it_sis(f1, n, res);
	for (int i = 0; i < n; i++)
	{
		cout << res[i] << endl;
	}*/
	/*int n = 4;
   double** mat1 = new double* [n];
   for (int i = 0; i < n; i++)
   {
	   mat1[i] = new double[n+1];
   }
   mat1[0][0] = -6;
   mat1[0][1] = -8;
   mat1[0][2] = -2;
   mat1[0][3] = -8;
   mat1[0][4] = -32;
   mat1[1][0] = 9;
   mat1[1][1] = 0;
   mat1[1][2] = 8;
   mat1[1][3] = 3;
   mat1[1][4] = 8;
   mat1[2][0] = 0;
   mat1[2][1] = -9;
   mat1[2][2] = -5;
   mat1[2][3] = 9;
   mat1[2][4] = -2;
   mat1[3][0] = -1;
   mat1[3][1] = 4;
   mat1[3][2] = -8;
   mat1[3][3] = -4;
   mat1[3][4] = -36;
   gayss(mat1, n, n + 1);
   for (int i = 0; i < n; i++)
   {
	   cout << mat1[i][n] << "  ";
   }*/
   /*int n = 4;
	   double** mat = new double* [n];
	   double**L = new double*[n];
	   double** U = new double* [n];
	   double* b = new double[n];
	   double* res = new double[n];
	   for (int i = 0; i < n; i++)
	   {
		   mat[i] = new double[n];
		   L[i] = new double[n];
		   U[i] = new double[n];
	   }
	   mat[0][0] = -6;
	   mat[0][1] = -8;
	   mat[0][2] = -2;
	   mat[0][3] = -8;
	   b[0] = -32;
	   mat[1][0] = 9;
	   mat[1][1] = 0;
	   mat[1][2] = 8;
	   mat[1][3] = 3;
	   b[1] = 8;
	   mat[2][0] = 0;
	   mat[2][1] = -9;
	   mat[2][2] = -5;
	   mat[2][3] = 9;
	   b[2] = -2;
	   mat[3][0] = -1;
	   mat[3][1] = 4;
	   mat[3][2] = -8;
	   mat[3][3] = -4;
	   b[3] = -36;
	   cout << LU_clay(mat, b, res, n) << endl;
	   for (int i = 0; i < n; i++)
	   {
		   cout << res[i] << "  ";
	   }
	   for (int i = 0; i < n; i++)
	   {
		   delete[]mat[i];
	   }
	   delete[]mat;*/
	   /*int n = 5;
		   double** mat = new double* [n];
		   double* res = new double[n];
		   for (int i = 0; i < n; i++)
		   {
			   mat[i] = new double[n+1];
			   for (int j = 0; j < n; j++)
			   {
				   mat[i][j] = 0;
			   }
		   }
		   mat[0][0] = 7;
		   mat[0][1] = -2;
		   mat[1][0] = -3;
		   mat[1][1] = -7;
		   mat[1][2] = 4;
		   mat[2][1] = -2;
		   mat[2][2] = 15;
		   mat[2][3] = 5;
		   mat[3][2] = -2;
		   mat[3][3] = -12;
		   mat[3][4] = -8;
		   mat[4][3] = -3;
		   mat[4][4] = -10;
		   mat[0][n] = 65;
		   mat[1][n] = 23;
		   mat[2][n] = 1;
		   mat[3][n] = -58;
		   mat[4][n] = -8;
		   res = progonka(mat, n);
		   for (int i = 0; i< n; i++)
		   {
			   cout << res[i] << endl;
		   }*/
		   //int n = 5;
			   //double** mat = new double* [n];
			   //double* res = new double[n];
			   //double* b = new double[n];
			   //for (int i = 0; i < n; i++)
			   //{
			   //	mat[i] = new double[n];
			   //	for (int j = 0; j < n; j++)
			   //	{
			   //		mat[i][j] = 0;
			   //	}
			   //}
			   //mat[0][0] = 7;
			   //mat[0][1] = -2;
			   //mat[1][0] = -3;
			   //mat[1][1] = -7;
			   //mat[1][2] = 4;
			   //mat[2][1] = -2;
			   //mat[2][2] = 15;
			   //mat[2][3] = 5;
			   //mat[3][2] = -2;
			   //mat[3][3] = -12;
			   //mat[3][4] = -8;
			   //mat[4][3] = -3;
			   //mat[4][4] = -10;
			   //b[0] = 65;
			   //b[1] = 23;
			   //b[2] = 1;
			   //b[3] = -58;
			   //b[4] = -8;
			   //res[0]=6.5;//7
			   //res[1]=-7;//-8
			   //res[2]=-2.5;//-3
			   //res[3]=5.5;//6
			   //res[4]=-1.5;//-1
			   //m_Seidel(mat, b, n, res);
			   //for (int i = 0; i < n; i++)
			   //{
			   //	cout << res[i] << endl;
			   //}
//int n = 5;
	//double** mat = new double* [n];
	//double* res = new double[n];
	//double* b = new double[n];
	//for (int i = 0; i < n; i++)
	//{
	//	mat[i] = new double[n];
	//	for (int j = 0; j < n; j++)
	//	{
	//		mat[i][j] = 0;
	//	}
	//}
	//mat[0][0] = 7;
	//mat[0][1] = -2;
	//mat[1][0] = -3;
	//mat[1][1] = -7;
	//mat[1][2] = 4;
	//mat[2][1] = -2;
	//mat[2][2] = 15;
	//mat[2][3] = 5;
	//mat[3][2] = -2;
	//mat[3][3] = -12;
	//mat[3][4] = -8;
	//mat[4][3] = -3;
	//mat[4][4] = -10;
	//b[0] = 65;
	//b[1] = 23;
	//b[2] = 1;
	//b[3] = -58;
	//b[4] = -8;
	//res[0] = 6.5;//7
	//res[1] = -7;//-8
	//res[2] = -2.5;//-3
	//res[3] = 5.5;//6
	//res[4] = -1.5;//-1
	//m_Jacobi(mat, b, n, res);
	//for (int i = 0; i < n; i++)
	//{
	//	cout << res[i] << endl;
	//}
/*int n = 4;
	double omega;
	cin>>omega;
	double** mat = new double* [n];
	double* res = new double[n];
	double* b = new double[n];
	for (int i = 0; i < n; i++)
	{
		mat[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			mat[i][j] = 0;
		}
	}
	mat[0][0] = 24;
	mat[0][1] = 9;
	mat[0][2] = -1;
	mat[0][3] = -5;
	mat[1][0] = -1;
	mat[1][1] = -14;
	mat[1][2] = 1;
	mat[1][3] = 9;
	mat[2][0] = -7;
	mat[2][1] = 5;
	mat[2][2] = -21;
	mat[3][0] = 1;
	mat[3][1] = 4;
	mat[3][2] = 8;
	mat[3][3] = -22;
	b[0] = -24;
	b[1] = 40;
	b[2] = -84;
	b[3] = -56;
	res[0] = 6.5;
	res[1] = -7;
	res[2] = -2.5;
	res[3] = 5.5;
	pvr(mat, b, n, res, Omega);
	for (int i = 0; i < n; i++)
	{
		cout << res[i] << endl;
	}
	*/
	//int n;
	//n = 3;
	//double** mat = new double* [n];
	//double** res = new double* [n];
	//double* lambda = new double[n];
	//for (int i = 0; i < n; i++)
	//{
	//	mat[i] = new double[n];
	//	res[i] = new double[n];
	//}
	//mat[0][0] = -4;
	//mat[0][1] = 1;
	//mat[0][2] = 7;
	//mat[1][0] = 1;
	//mat[1][1] = 8;
	//mat[1][2] = -5;
	//mat[2][0] = 7;
	//mat[2][1] = -5;
	//mat[2][2] = 1;
	//m_rot(mat, n, res, lambda);
	//for (int i = 0; i < n; i++)
	//{
	//	for (int j = 0; j < n; j++)
	//	{
	//		cout << res[j][i] << "   ";
	//	}
	//	cout << lambda[i];
	//	cout << endl;
	//}
////int n=3;
	////double** mat = new double* [n];
	////for (int i = 0; i < n; i++)
	////{
	////	mat[i] = new double[n];
	////}
	////mat[0][0] = -4;
	////mat[0][1] = 1;
	////mat[0][2] = 7;
	////mat[1][0] = 1;
	////mat[1][1] = 8;
	////mat[1][2] = -5;
	////mat[2][0] = 7;
	////mat[2][1] = -5;
	////mat[2][2] = 1;
	////cout<<max_lambda(mat, n);
//int n;
//	n = 3;
//	double** mat = new double* [n];
//	double** res = new double* [n];
//	double* lambda = new double[n];
//	for (int i = 0; i < n; i++)
//	{
//		mat[i] = new double[n];
//		res[i] = new double[n];
//	}
//	mat[0][0] = -4;
//	mat[0][1] = 1;
//	mat[0][2] = 7;
//	mat[1][0] = 1;
//	mat[1][1] = 8;
//	mat[1][2] = -5;
//	mat[2][0] = 7;
//	mat[2][1] = -5;
//	mat[2][2] = 1;
//	m_QR(mat, n, lambda, res);
//	for (int i = 0; i < n; i++)
//	{
//		cout << lambda[i];
//		cout << endl;
//	}
//double h1, h2,a,b;
//h1 = 1;
//h2 = 0.5;
//a = 1;
//b = 5;
//cout << S_pr(f4, a, b, h1) << endl;
//cout << S_pr(f4, a, b, h2) << endl;
//cout << S_trap(f4, a, b, h1) << endl;
//cout << S_trap(f4, a, b, h2) << endl;
//cout << S_Simpson(f4, a, b, h1) << endl;
//cout << S_Simpson(f4, a, b, h2) << endl;
//cout << Romberg(f4, a, b);
    return 0;
}