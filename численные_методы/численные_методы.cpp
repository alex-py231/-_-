using namespace std;
#include <iostream>
#include "функции.h"
#define Omega 1.
int main()
{
    setlocale(LC_ALL, "Russian");
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

	//int n = 4;
	//double omega = 0.99;
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
	//mat[0][0] = 24;
	//mat[0][1] = 9;
	//mat[0][2] = -1;
	//mat[0][3] = -5;
	//mat[1][0] = -1;
	//mat[1][1] = -14;
	//mat[1][2] = 1;
	//mat[1][3] = 9;
	//mat[2][0] = -7;
	//mat[2][1] = 5;
	//mat[2][2] = -21;
	//mat[3][0] = 1;
	//mat[3][1] = 4;
	//mat[3][2] = 8;
	//mat[3][3] = -22;
	//b[0] = -24;
	//b[1] = 40;
	//b[2] = -84;
	//b[3] = -56;
	//res[0] = 6.5;
	//res[1] = -7;
	//res[2] = -2.5;
	//res[3] = 5.5;
	//pvr(mat, b, n, res, Omega);
	//for (int i = 0; i < n; i++)
	//{
	//	cout << res[i] << endl;
	//}
	/*double x=0.8;
	int n=4;
	double* X = new double [n];
	double* Y = new double[n];
	X[0] = 0.1;
	X[1] = 0.5;
	X[2] = 0.9;
	X[3] = 1.3;
	Y[0] = 0.1;
	Y[1] = 0.5;
	Y[2] = 1.1;
	Y[3] = 1.3;
	cout<<inter_laplas(X, Y, n, x);*/
	 
	 
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
	/*int n=3;
	double** mat = new double* [n];
	for (int i = 0; i < n; i++)
	{
		mat[i] = new double[n];
	}
	mat[0][0] = -4;
	mat[0][1] = 1;
	mat[0][2] = 7;
	mat[1][0] = 1;
	mat[1][1] = 8;
	mat[1][2] = -5;
	mat[2][0] = 7;
	mat[2][1] = -5;
	mat[2][2] = 1;
	double* lambda = new double[n];
	lambda[0] = 11.0788;
	lambda[1] = 3.6607;
	lambda[2] = -9.73954;
	double** mat1 = new double* [n];
	for (int i = 0; i < n; i++)
	{
		mat1[i] = new double[n+1];
	}
	for (int k = 0; k < n; k++)
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (i == j)
				{
					mat1[i][i] = mat[i][i] - lambda[k];
				}
				else 
				{
					mat1[i][j] = mat[i][j];
				}
			}
			mat1[i][n] = 0;
		}
		gayss(mat1, n, n + 1);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				cout << mat1[i][j]<<"   ";
			}
			cout << endl;
		}
		cout << endl;
	}*/
	int n;
	n=3;
	double** mat = new double*[n];
	double** res = new double* [n];
	for (int i = 0; i < n; i++)
	{
		mat[i] = new double[n];
		res[i] = new double[n];
	}
	mat[0][0] = -4;
	mat[0][1] = 1;
	mat[0][2] = 7;
	mat[1][0] = 1;
	mat[1][1] = 8;
	mat[1][2] = -5;
	mat[2][0] = 7;
	mat[2][1] = -5;
	mat[2][2] = 1;
	m_rot(mat, n, res);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << res[j][i] << "   ";
		}
		cout << endl;
	}
    return 0;
}