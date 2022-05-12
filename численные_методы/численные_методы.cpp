using namespace std;
#include <iostream>
#include "функции.h"
#define Omega 1.
double f1(double x)
{
	return sqrt(x) / (4 + 3 * x);
}
int main()
{
    setlocale(LC_ALL, "Russian");
	//int n;
	//n=3;
	//double** mat = new double* [n];
	//double** mat1 = new double* [n];
	//double* lambda = new double[n];
	//double** sob_vec = new double* [n];
	//for (int i = 0; i < n; i++)
	//{
	//	mat[i] = new double[n];
	//	mat1[i] = new double[n];
	//	sob_vec[i] = new double[n];
	//}
	//mat[0][0] = 7;
	//mat[0][1] = 4;
	//mat[0][2] = 2;
	//mat[1][0] = -8;
	//mat[1][1] = -6;
	//mat[1][2] = 5;
	//mat[2][0] = 9;
	//mat[2][1] = 0;
	//mat[2][2] = -8;
	//m_rot(mat, n, sob_vec, lambda);       
	//for (int i = 0; i < n; i++)
	//{
	//	for (int j = 0; j < n; j++)
	//	{
	//		cout << sob_vec[i][j] << "   ";
	//	}
	//	cout << lambda[i];
	//	cout << endl;
	//}
	///*double a, b, h1, h2;
	//a = 1;
	//b = 5;
	//h1 = 1;
	//h2 = 0.5;
	//cout << h1 << ":" << S_pr(f1, a, b, h1)<<endl;
	//cout << h2 << ":" << S_pr(f1, a, b, h2) << endl;
	//cout << h1 << ":" << S_trap(f1, a, b, h1) << endl;
	//cout << h2 << ":" << S_trap(f1, a, b, h2) << endl;
	//cout << h1 << ":" << S_Simpson(f1, a, b, h1) << endl;
	//cout << h2 << ":" << S_Simpson(f1, a, b, h2) << endl;
	//cout << Romberg(f1, a, b);*/
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

    return 0;
}