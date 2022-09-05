using namespace std;
#include <iostream>
#include "функции.h"
int sign(int i, int j)
{
	if(i == j)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
double f_i(double x, double y,double z,double psi,double ksi,int i)
{
	double res = 0.1 * sinh(x) + pow((x - y), 3) + pow((x - z), 3) + 0.11 * (x - psi) + 0.11 * (x - ksi) + sign(i, 0) * (x - 11.) + sign(i, 23) * (x + 11);
	return res;
}

typedef double (*F)(double x, double y, double z, double psi, double ksi, int i);
typedef double (*dF)(F f, double x, double y, double z, double psi, double ksi, int i);
double df_x(F f, double x, double y, double z, double psi, double ksi, int i)
{
	double e = 0.000001;
	return ((f(x + e, y, z, psi, ksi, i) - f(x, y, z, psi, ksi, i)) / e);
}
double df_y(F f, double x, double y, double z, double psi, double ksi, int i)
{
	double e = 0.000001;
	return (f(x , y+e, z, psi, ksi, i) - f(x, y, z, psi, ksi, i)) / e;
}
double df_z(F f, double x, double y, double z, double psi, double ksi, int i)
{
	double e = 0.000001;
	return (f(x , y, z+e, psi, ksi, i) - f(x, y, z, psi, ksi, i)) / e;
}
double df_psi(F f, double x, double y, double z, double psi, double ksi, int i)
{
	double e = 0.000001;
	return (f(x , y, z, psi+e, ksi, i) - f(x, y, z, psi, ksi, i)) / e;
}
double df_ksi(F f, double x, double y, double z, double psi, double ksi, int i)
{
	double e = 0.000001;
	return (f(x , y, z, psi, ksi+e, i) - f(x, y, z, psi, ksi, i)) / e;
}
void newton_sis_T(F* f_, dF*df, int dim, double* res)
{
	double** mat = new double* [dim];
	double eps = 0.000001;
	double xmax ;
	double* dx = new double[dim];
	double x; double y; double z; double psi; double ksi;
	for (int i = 0; i < dim; i++)
	{
		mat[i] = new double[dim + 1];
	}
	do {
		xmax = 0.000001;
		int j, k, l, n;
		for (int i = 0; i < dim; i++)
		{
			j = i - 1;
			if (j < 0)
			{
				j = 50 - i - 1;
			}
			k = i + 1;
			if (k > 49)
			{
				k = 50 - i + 1;
			}
			l = i - 13;
			if (l < 0)
			{
				l = 50 + (i - 13);
			}
			n = i + 13;
			if (n > 49)
			{
				n = i + 13 - 50;
			}
			x = res[i];
			y = res[j];
			z = res[k];
			psi = res[l];
			ksi = res[n];
			for (int p = 0;p < dim;p++)
			{
				if (p != i && p != j && p != k && p != l && p != n)
				{
					mat[i][p] = 0;
				}
				else
				{
					if (p == i)
					{
						mat[i][p] = df[0](f_i, x, y, z, psi, ksi, i);
					}
					else if (p == j)
					{
						mat[i][j] = df[1](f_i, x, y, z, psi, ksi, i);
					}
					else if (p == k)
					{
						mat[i][p] = df[2](f_i, x, y, z, psi, ksi, i);
					}
					else if (p == l)
					{
						mat[i][p] = df[3](f_i, x, y, z, psi, ksi, i);
					}
					else
					{
						mat[i][p] = df[4](f_i, x, y, z, psi, ksi, i);
					}
				}
			}
			mat[i][dim] = -(f_i(x, y, z, psi, ksi, i));
		}
		gayss(mat, dim, dim + 1);
		for (int i = 0; i < dim; i++)
		{
			dx[i] = mat[i][dim];
		}
		for (int i = 0; i < dim; i++)
		{
			res[i] += dx[i];
		}
		for (int i = 0; i < dim; i++)
		{
			xmax = (fabs(dx[i]) > xmax) ? fabs(dx[i]) : xmax;
		}
	} while (fabs(xmax) > eps);
}
int main()
{
    setlocale(LC_ALL, "Russian");
	int n = 50;
	F* f1 = new F[n];
	dF* df = new dF[5];
	double* res = new double[n];
	for (int i = 0;i < n;i++)
	{
		res[i] = 0;
	}
	for (int i = 0;i < n;i++)
	{
		f1[i] = f_i;
	}
	df[0] = df_x;
	df[1] = df_y;
	df[2] = df_z;
	df[3] = df_psi;
	df[4] = df_ksi;
	newton_sis_T(f1, df, n, res);
	for (int i = 0; i < n; i++)
	{
		cout << res[i] << endl;
	}

	/*int j, k, l, p;
	double x, y, z, ksi, psi;
	for (int i = 0; i < n; i++)
	{
		j = i - 1;
		if (j < 0)
		{
			j = 50 - i - 1;
		}
		k = i + 1;
		if (k > 49)
		{
			k = 50 - i + 1;
		}
		l = i - 13;
		if (l < 0)
		{
			l = 50 + (i - 13);
		}
		p = i + 13;
		if (p > 49)
		{
			p = i + 13 - 50;
		}
		x = res[i];
		y = res[j];
		z = res[k];
		psi = res[l];
		ksi = res[p];
		cout << (f_i(x, y, z, psi, ksi, i)) << endl;
	}*/     // проверка результатов 

    return 0;
}