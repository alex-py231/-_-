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





	int n = 5;
	double omega = 0.99;
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
	b[0] = 65;
	b[1] = 23;
	b[2] = 1;
	b[3] = -58;
	b[4] = -8;
	res[0] = 6.5;//7
	res[1] = -7;//-8
	res[2] = -2.5;//-3
	res[3] = 5.5;//6
	res[4] = -1.5;//-1
	pvr(mat, b, n, res, Omega);
	for (int i = 0; i < n; i++)
	{
		cout << res[i] << endl;
	}
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

    return 0;
}