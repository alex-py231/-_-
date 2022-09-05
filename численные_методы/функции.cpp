#include "функции.h"
double** traspose(double** mat, int n, int m)
{
	double** arr = new double* [n];
	for (int i = 0; i < n; i++)
	{
		arr[i] = new double[m];
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			arr[i][j] = mat[j][i];
		}
	}
	return arr;
}
double helpfunk(double** arr, int n, int a)
{
	for (int i = 0; i < n; i++)
	{
		if (arr[i][a] == 0);
		else
		{
			return i;
			break;
		}
	}
}
void obr(double** mat, int n)
{
	double** arr = new double* [n];
	for (int i = 0; i < n; i++)
	{
		arr[i] = new double[2 * n];
		for (int j = 0; j < n; j++)
		{
			arr[i][j] = mat[i][j];
		}
		for (int j = n; j < 2 * n; j++)
		{
			if (i == j - n)
			{
				arr[i][j] = 1;
			}
			else
			{
				arr[i][j] = 0;
			}
		}
	}
	for (int a = 0; a < n; a++)
	{
		{
			if (arr[a][a] == 0)
			{
				int i = helpfunk(arr, n, a);
				for (int j = 0; j < n; j++)
				{
					arr[a][j] = arr[a][j] + arr[i][j];
				}
			}
			double c = arr[a][a];
			for (int j = 0; j < 2 * n; j++)
			{
				arr[a][j] = (double)arr[a][j] / c;
			}
			for (int i = 0; i < n; i++)
			{
				double d = arr[i][a];
				if (i == a);
				else
				{
					if (arr[i][a] == 0);
					if (arr[i][a] != 0)
					{
						for (int j = 0; j < 2 * n; j++)
						{
							arr[i][j] = (double)arr[i][j] - d * arr[a][j];
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = n; j < 2 * n; j++)
		{
			mat[i][j - n] = arr[i][j];
		}
	}
}
double L(double* X, int n, double x, int i)//вспомогательная функция к функции Лапласа
{
	double res = 1;
	for (int j = 0; j < n; j++)
	{
		if (j != i)
		{
			res *= (x - X[j]);
		}
		else
		{
			res *= 1;
		}
	}
	return res;
}
double l(double* X, int n, int i)//Вспомогательная функция к функции Лапласа
{
	double res = 1;
	for (int j = 0; j < n; j++)
	{
		if (j != i)
		{
			res *= (X[i] - X[j]);
		}
		else
		{
			res *= 1;
		}
	}
	return res;
}
double inter_laplas(double* X, double* Y, int n, double x)//функция интерполяции методом Лапласа(2 лекция) (вспомогательные к ней L,l)
{//идея задать значения независимой переменой и значения функции от этой переменой ,а промежуточные посчитает функция
	double P = 0;
	for (int i = 0; i < n; i++)
	{
		P += Y[i] * L(X, n, x, i) / l(X, n, i);
	}
	return P;
}
double P_(double* X, double* Y, int start, int stop)
{
	if (stop - start == 1)
	{
		return ((Y[stop] - Y[start]) / (X[stop] - X[start]));
	}
	else if (stop - start > 1)
	{
		return ((P_(X, Y, start + 1, stop) - P_(X, Y, start, stop - 1)) / (X[stop] - X[start]));
	}
	else if (stop - start == 0)
	{
		return 1;
	}
}
double P(double* X, double* Y, int n)
{
	double pr;
	double res = 0;
	for (int i = 0; i <= n; i++)
	{
		pr = 1;
		for (int j = 0; j <= n; j++)
		{
			if (j != i)
			{
				pr *= (X[i] - X[j]);
			}
		}
		res += Y[i] / pr;
	}
	return res;
}
double Pr(double* X, int k, double x)
{
	double res = 1;
	for (int i = 0; i < k; i++)
	{
		res *= (x - X[i]);
	}
	return res;
}
double inter_Newton(double* X, double* Y, int n, double x)
{//идея интерполировать функцию конечными разностями(P(),P_())
	double res = Y[0];
	for (int i = 1; i < n; i++)
	{
		res += (P(X, Y, i) * Pr(X, i, x));
	}
	return res;
}
double error_rate(double y, double Y)//Y - точное значение,y- приближёное значение
{
	return abs(Y - y) / abs(Y) * 100;
}
void coef(double* Y, double h, int N, double* b_coef, double* c, double* d)
{
	c[0] = 0;
	double* a = new double[N];
	double* b = new double[N];
	double F, A, B, C, z;
	a[0] = 0;
	b[0] = 0;
	for (int i = 1; i < N - 1; i++)
	{
		A = h;
		C = 2. * 2 * h;
		B = h;
		F = 6. * ((Y[i + 1] - Y[i]) / h - (Y[i] - Y[i - 1]) / h);
		z = (A * a[i - 1] + C);
		a[i] = -B / z;
		b[i] = (F - A * b[i - 1]) / z;
	}
	c[N - 1] = (F - A * b[N - 2]) / (C + A * a[N - 2]);
	for (int i = N - 2; i > 0; i--)
	{
		c[i] = a[i] * c[i + 1] - b[i];
	}
	delete[]a;
	delete[]b;
	for (int i = N - 1; i > 0; --i)
	{
		d[i] = (c[i] - c[i - 1]) / h;
		b_coef[i] = h * (2. * c[i] - c[i - 1]) / 6. + (Y[i] - Y[i - 1]) / h;
	}

}
int search_max(double** arr, int n, int j)
{
	double max = arr[0][j];
	int i_max = 0;
	for (int i = 0; i < n; i++)
	{
		if (max < abs(arr[i][j]))
		{
			max = arr[i][j];
			i_max = i;
		}
	}
	return i_max;
}
void swap_row(double** arr, int i, int j, int n)
{
	double* a = new double[n];
	for (int k = 0; k < n; k++)
	{
		a[k] = arr[i][k];
		arr[i][k] = arr[j][k];
		arr[j][k] = a[k];
	}
	delete[]a;
}
void epsilon(double** arr, int n, int m, double eps)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (abs(arr[i][j]) < eps)
			{
				arr[i][j] = 0;
			}
		}
	}
}
void gayss(double** arr, int n, int m)
{
	int a = 0;
	int i;
	double eps = 0.00001;
	do {
		epsilon(arr, n, m, eps);
		i = search_max(arr, n, a);
		if (i == a);
		else
		{
			swap_row(arr, a, i, m);
		}
		double c = arr[a][a];//константа 1
		for (int j = 0; j < m; j++)
		{
			arr[a][j] = (double)arr[a][j] / c;
		}
		for (int i = 0; i < n; i++)
		{
			double c1 = arr[i][a];//констата 2
			if (i == a);
			else
			{
				if (abs(arr[i][a]) < eps);
				if (abs(arr[i][a]) > eps)
				{
					for (int j = 0; j < m; j++)
					{
						arr[i][j] = (double)arr[i][j] - c1 * arr[a][j];
					}
				}
			}
		}
		a++;
	} while (a < n);
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
}
double power(int n, double x)
{
	double res = 1;
	for (int i = 0; i < n; i++)
	{
		res *= x;
	}
	return res;
}
void MNK(double* X, double* Y, double* res, int n, int m)
{//аппраксимация полиномом любой степени методом наимньших квадратов.
	double** mat = new double* [n];
	for (int i = 0; i < n; i++)
	{
		mat[i] = new double[n + 1];
		for (int j = 0; j < n + 1; j++)
		{
			mat[i][j] = 0.;
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < m; k++)
			{
				mat[i][j] += power(j + i, X[k]);
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int k = 0; k < m; k++)
		{
			mat[i][n] += Y[k] * power(i, X[k]);
		}
	}
	gayss(mat, n, n + 1);
	for (int i = 0; i < n; i++)
	{
		res[i] = mat[i][n];
	}
}
double cubic_spline(double* X, double* Y, int n, double x)
{//какая-то шляпа ,код явно нужно переписать вся пробелема вроде бы в функции coef
//когда будешь переписывать используй нормальный метод прогонки,он должен быть ниже
	double h = X[1] - X[0];
	int s;
	double* b = new double[n];
	double* c = new double[n];
	double* d = new double[n];
	coef(Y, h, n, b, c, d);
	if (x <= X[0])
	{
		s = 1;
	}
	else if (x >= X[n - 1])
		s = n - 1;
	else
	{
		int i = 0, j = n - 1;
		while (i + 1 < j)
		{
			int k = i + (j - i) / 2;
			if (x <= X[k])
				j = k;
			else
				i = k;
		}
		s = j - 1;
	}
	double dx = (x - X[s]);
	return (Y[s] + (b[s] + (c[s] / 2. + d[s] * dx / 6.) * dx) * dx) / 2;
}
double* progonka(double** mat, int n)//где размер матрицы n на n+1
{// !этот метод применяется только к 3-х диоганальым матрицам!
	double* a = new double[n];
	double* b = new double[n];
	double* res = new double[n];
	double F, A, B, C, z;
	A = 0;
	a[0] = 0;
	b[0] = 0;
	for (int i = 0; i < n - 1; i++)
	{
		if (i != 0)
		{
			A = mat[i][i - 1];
		}
		B = mat[i][i];
		C = mat[i][i + 1];
		F = mat[i][n];
		a[i + 1] = -C / (A * a[i] + B);
		b[i + 1] = (F - A * b[i]) / (A * a[i] + B);
	}
	res[n - 1] = (mat[n - 1][n] - mat[n - 1][n - 2] * b[n - 1]) / (mat[n - 1][n - 1] + mat[n - 1][n - 2] * a[n - 1]);
	for (int i = n - 2; i >= 0; i--)
	{
		res[i] = a[i + 1] * res[i + 1] + b[i + 1];
	}
	delete[]a;
	delete[]b;
	return res;

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
}
double m_hord(double (*f)(double), double a, double b)//метод хорд- приближёное решение уравнения,где
//a,b точки отрезка где ищем корень,f -функция уравнения.eps - точность.
{
	double eps = 0.0001;
	do {
		a = b - (b - a) * f(b) / (f(b) - f(a));
		b = a - (a - b) * f(a) / (f(a) - f(b));
	} while (abs(b - a) > eps);
	return b;
}
double m_newtona(double (*fx)(double), double (*dfx)(double), double x0) {
	double x1 = x0 - fx(x0) / dfx(x0); // первое приближение
	double eps = 0.0001;
	while (abs(x1 - x0) > eps) { // пока не достигнута точность 0.0001
		x0 = x1;
		x1 = x0 - fx(x0) / dfx(x0); // последующие приближения
	}
	return x1;
}
double m_pros_it(double (*f)(double), double a, double b)
{
	double x_k, x_k1 = 0;
	double eps = 0.0001;
	x_k = b;
	while (abs(x_k - x_k1) > eps) {
		x_k1 = x_k;
		x_k = f(x_k1);
	};
	return x_k;
}
bool LU(double** mat, int n, int m, double** L, double** U)
{
	double eps = 0.00001;
	double** tmp = new double* [n];
	for (int i = 0; i < n; i++)
	{
		tmp[i] = new double[2 * m];
		for (int j = 0; j < m; j++)
		{
			tmp[i][j] = mat[i][j];
		}
		for (int j = m; j < 2 * m; j++)
		{
			if (i == j - m)
			{
				tmp[i][j] = 1;
			}
			else
			{
				tmp[i][j] = 0;
			}
		}
	}
	double tmp_;
	for (int k = 0; k < n - 1; k++)
	{
		epsilon(tmp, n, 2 * n, eps);
		if (abs(tmp[k][k]) < eps)
		{
			return false;
		}
		else
		{
			for (int i = k + 1; i < n; i++)
			{
				tmp_ = tmp[i][k] / tmp[k][k];
				for (int j = 0; j < 2 * m; j++)
				{
					tmp[i][j] -= tmp_ * tmp[k][j];
				}
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			U[i][j] = tmp[i][j];
		}
		for (int j = m; j < 2 * m; j++)
		{
			L[i][j - m] = tmp[i][j];
		}
	}
	obr(L, n);
	return true;
}
double** pr_mat(double** a, double** b, int n, int m)
{
	double** res = new double* [n];
	for (int j = 0; j < n; j++)
	{
		res[j] = new double[m];
		for (int i = 0; i < m; i++)
		{
			res[j][i] = 0;
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			for (int k = 0; k < m; k++)
			{
				res[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	return res;
}
void pr_mat_2(double** a, double** b, int n, int m, double** res1)
{
	double** res = new double* [n];
	for (int j = 0; j < n; j++)
	{
		res[j] = new double[m];
		for (int i = 0; i < m; i++)
		{
			res[j][i] = 0;
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			for (int k = 0; k < m; k++)
			{
				res[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			res1[i][j] = res[i][j];
		}
	}
}
bool LU_clay(double** mat, double* b, double* res, int n)// решение СЛАУ происходит так ,в начале матрица коэффициентов разлагается на верхнетреугольную и нижнетреугольную,а дальше вектор столбец b доможается на L-1,u-1
{
	double** L = new double* [n];
	double** U = new double* [n];
	double* tmp = new double[n];
	double tmp_ = 0;
	for (int i = 0; i < n; i++)
	{
		L[i] = new double[n];
		U[i] = new double[n];
	}
	if (LU(mat, n, n, L, U))
	{
		for (int i = 0; i < n; i++)//обратный ход метода Гаусса
		{
			tmp_ = 0;
			for (int j = 0; j < i; j++)
			{
				tmp_ += tmp[j] * L[i][j];
			}
			tmp[i] = (b[i] - tmp_) / L[i][i];
		}
		for (int i = n - 1; i >= 0; i--)//обратный ход метода Гаусса
		{
			tmp_ = 0;
			for (int j = i + 1; j < n; j++)
			{
				tmp_ += res[j] * U[i][j];
			}
			res[i] = (tmp[i] - tmp_) / U[i][i];
		}
		return true;
	}
	else
	{
		return false;
	}
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
}
//void pr_it_sis(f* f_, int n, double* res)//res- начальное приблежение,в нём вернётся и результат, *f_ - массив функций
//{// по понятным причинам n -  размер массива
//	double eps = 0.0001;
//	double* res1 = new double[n];
//	double sum = 0;
//	double dx;
//	do {
//		dx = 0.00001;
//		for (int i = 0; i < n; i++)
//		{
//			res1[i] = f_[i](res);
//		}
//		for (int i = 0; i < n; i++)
//		{
//			dx = (abs((res[i] - res1[i])) > dx) ? abs(res[i] - res1[i]) : dx;
//		}
//		for (int i = 0; i < n; i++)
//		{
//			res[i] = res1[i];
//		}
//	} while (dx > eps );
//	delete[]res1;
//	//double f1(double* x)
//		//{
//		//    return sqrt(2 * x[0] * x[0] - x[0] - 1);
//		//}
//		//double f2(double *x)
//		//{
//		//    return atan(x[1]);
//		//}
//		/* f* f_ = new f[n];
//	f_[0] = f1;
//	f_[1] = f2;
//	double* res = new double[n];
//	res[0] = 0.6;
//	res[1] = 0.6;
//	pr_it_sis(f_, n, res);
//	for (int i = 0; i < n; i++)
//	{
//		cout << res[i] << "   ";
//	}*/
//}
//void newton_sis(f* f_, f** df, int n, double* res)// по параметрам,та же история,что и в прошлый раз
//{// df двумерный массив дифференцированных функций,фактически матрица Якоби
//	double** mat = new double* [n];
//	double eps = 0.0001;
//	double xmax;
//	double *dx = new double[n];
//	for (int i = 0; i < n; i++)
//	{
//		mat[i] = new double[n+1];
//	}
//	do {
//		xmax = 0.00001;
//		for (int i = 0; i < n; i++)
//		{
//			for (int j = 0; j < n; j++)
//			{
//				mat[i][j] = df[i][j](res);
//			}
//			mat[i][n] = -(f_[i](res));
//		}
//		gayss(mat, n, n + 1);
//		for (int i = 0; i < n; i++)
//		{
//			dx[i] = mat[i][n];
//		}
//		for (int i = 0; i < n; i++)
//		{
//			res[i] += dx[i];
//		}
//		for (int i = 0; i < n; i++)
//		{
//			xmax = (fabs(dx[i]) > xmax) ? fabs(dx[i]) : xmax;
//		}
//	} while (fabs(xmax) > eps);
//	//
////double f1(double* x)
////{
////    return x[0] * x[0] + 2 * log10(x[1]) - 1;
////}
////double f2(double *x)
////{
////    return x[0] * x[0] - 3 * x[0] * x[1] + 3;
////}
////double df1(double* x)
////{
////    return 2 * x[0];
////}
////double df2(double* x)
////{
////    return 2 / (x[1] * log(10));
////}
////double df3(double* x)
////{
////    return 2 * x[0] - 3 * x[1];
////}
////double df4(double* x)
////{
////    return -3 * x[0];
////}
//	/* f** df = new f * [n];
//	for (int i = 0; i < n; i++)
//	{
//		df[i] = new f[n];
//	}
//	df[0][0] = df1;
//	df[0][1] = df2;
//	df[1][0] = df3;
//	df[1][1] = df4;
//	f* F = new f[n];
//	F[0] = f1;
//	F[1] = f2;
//	double* res = new double[n];
//	res[0] = 0.3;
//	res[1] = 3;
//	newton_sis(F, df, n, res);
//	for (int i = 0; i < n; i++)
//	{
//		cout << res[i] << "   ";
//	}
//	res[0] = 1;
//	res[1] = 3;
//	newton_sis(F, df, n, res);
//	for (int i = 0; i < n; i++)
//	{
//		cout << res[i] << "   ";
//	}*/
//}
void m_Seidel(double** mat, double* b, int n, double* res)//mat- матрица коэффициентов,b- вектор столбец свободных членов,res - первое приближение,в нём же вернётся и результат
{
	double* res1 = new double[n];
	double** C = new double* [n];
	double* d = new double[n];
	for (int i = 0; i < n; i++)
	{
		C[i] = new double[n];
		res1[i] = res[i];
	}
	for (int i = 0; i < n; i++)
	{
		d[i] = b[i] / mat[i][i];
		for (int j = 0; j < n; j++)
		{
			if (i == j) C[i][j] = 0;
			else
			{
				C[i][j] = -1. * mat[i][j] / mat[i][i];
			}
		}
	}
	double eps = 0.00001;
	double dx;
	do {
		dx = 0.000001;
		for (int i = 0; i < n; i++)
		{
			res[i] = 0;
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < i; j++)
			{
				res[i] += C[i][j] * res[j];
			}
			for (int j = i; j < n; j++)
			{
				res[i] += C[i][j] * res1[j];
			}
			res[i] += d[i];
		}
		for (int i = 0; i < n; i++)
		{
			dx = (fabs(res[i] - res1[i]) > dx) ? fabs(res[i] - res1[i]) : dx;
		}
		for (int i = 0; i < n; i++)
		{
			res1[i] = res[i];
		}
	} while (dx > eps);
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
}
void m_Jacobi(double** mat, double* b, int n, double* res)
{
	double* res1 = new double[n];
	for (int i = 0; i < n; i++)
	{
		res1[i] = res[i];
	}
	double eps = 0.00001;
	double dx;
	double sum;
	do {
		dx = 0.000001;
		for (int i = 0; i < n; i++)
		{
			res[i] = 0;
		}
		for (int i = 0; i < n; i++)
		{
			sum = 0;
			for (int j = 0; j < n; j++)
			{
				if (i != j)
				{
					sum += mat[i][j] * res1[j];
				}
				else continue;
			}
			res[i] = (b[i] - sum) / mat[i][i];
		}
		for (int i = 0; i < n; i++)
		{
			dx = (fabs(res[i] - res1[i]) > dx) ? fabs(res[i] - res1[i]) : dx;
		}
		for (int i = 0; i < n; i++)
		{
			res1[i] = res[i];
		}
	} while (dx > eps);
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
}
void pvr(double** mat, double* b, int n, double* res, double omega)
{
	double* res1 = new double[n];
	double** C = new double* [n];
	double* d = new double[n];
	for (int i = 0; i < n; i++)
	{
		C[i] = new double[n];
		res1[i] = res[i];
	}
	for (int i = 0; i < n; i++)
	{
		d[i] = b[i] / mat[i][i];
		for (int j = 0; j < n; j++)
		{
			if (i == j) C[i][j] = 0;
			else
			{
				C[i][j] = -mat[i][j] / mat[i][i];
			}
		}
	}
	double eps = 0.00001;
	double dx;
	do {
		dx = 0.000001;
		for (int i = 0; i < n; i++)
		{
			res[i] = 0;
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < i; j++)
			{
				res[i] += omega * C[i][j] * res[j];
			}
			for (int j = i; j < n; j++)
			{
				res[i] += omega * C[i][j] * res1[j];
			}
			res[i] += omega * d[i];
		}
		for (int i = 0; i < n; i++)
		{
			dx = (fabs(res[i] - res1[i]) > dx) ? fabs(res[i] - res1[i]) : dx;
		}
		for (int i = 0; i < n; i++)
		{
			res1[i] = res[i];
		}

	} while (dx > eps);

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
	// 
}
void m_rot(double** mat1, int n, double** res1, double* lambda)
{
	double** mat = new double* [n];
	double** res = new double* [n];
	double** H = new double* [n];
	double eps = 0.0001;
	for (int i = 0; i < n; i++)
	{
		mat[i] = new double[n];
		H[i] = new double[n];
		res[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			mat[i][j] = mat1[i][j];
			if (i == j)
			{
				res[i][j] = 1.;
			}
			else
			{
				res[i][j] = 0;
			}
		}
	}
	double max = 1;
	int i_max = 0, j_max = 0;
	double phi;
	while (max > eps) {
		max = 0;
		for (int i = 0; i < n; i++)
		{
			for (int j = i + 1; j < n; j++)
			{
				if (max < abs(mat[i][j]))
				{
					max = abs(mat[i][j]);
					i_max = i;
					j_max = j;
				}
			}
		}
		phi = atan(2 * mat[i_max][j_max] / (mat[i_max][i_max] - mat[j_max][j_max])) / 2;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (i == j)
				{
					H[i][j] = 1;
				}
				else
				{
					H[i][j] = 0;
				}
			}
		}
		H[i_max][i_max] = cos(phi);
		H[i_max][j_max] = -sin(phi);
		H[j_max][i_max] = sin(phi);
		H[j_max][j_max] = cos(phi);
		mat = pr_mat(traspose(H, n, n), mat, n, n);
		mat = pr_mat(mat, H, n, n);
		res = pr_mat(res, H, n, n);
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			res1[i][j] = res[i][j];
		}
		lambda[i] = mat[i][i];
	}
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
}
double S_pr(double(*f)(double), double a, double b, double h)
{
	int n = (b - a) / h;
	double sum = 0;
	for (double i = a; i < b; i += h)
	{
		sum += f((2 * i + h) / 2) * h;
	}
	return sum;
}
double S_trap(double(*f)(double), double a, double b, double h)
{
	int n = (b - a) / h;
	double sum = 0;
	for (double i = a + h; i < b; i += h)
	{
		sum += (f(i - h) + f(i)) / 2 * h;
	}
	return sum;
}
double S_Simpson(double(*f)(double), double a, double b, double h)
{
	int n = (b - a) / h;
	double sum = 0;
	for (double i = a + h; i < b; i += h)
	{
		sum += (f(i - h) + 4 * f(((i - h) + i) / 2) + f(i)) / 6 * h;
	}
	return sum;
}
double scal_pr(double* A, double* B, int n)
{
	double res = 0;
	for (int i = 0; i < n; i++)
	{
		res += A[i] * B[i];
	}
	return res;
}
void ort_gramm_shimidt(double** mat, double** res, int n)
{
	double* b = new double[n];
	double** mat1 = new double* [n];
	for (int i = 0; i < n; i++)
	{
		mat1[i] = new double[n];
	}
	double* c = new double[n];
	double sum = 0;
	double* m = new double[n];
	for (int j = 0; j < n; j++)
	{
		for (int k = 0; k < n; k++)
		{
			b[k] = mat[k][j];
		}
		for (int k = 0; k < j; k++)
		{
			for (int i = 0; i < n; i++)
			{
				m[i] = mat1[i][k];
			}
			c[k] = scal_pr(b, m, n) / scal_pr(m, m, n);
		}
		for (int k = 0; k < n; k++)
		{
			sum = 0;
			for (int m = 0; m < j; m++)
			{
				sum += c[m] * mat1[k][m];
			}
			b[k] -= sum;
			sum = 0;
		}
		for (int k = 0; k < n; k++)
		{
			mat1[k][j] = b[k];
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			res[i][j] = mat1[i][j];
		}
	}
}
void norm(double** mat, double** res, int n)
{
	double* b = new double[n];
	double mod = 0;
	for (int i = 0; i < n; i++)
	{
		mod = 0;
		for (int j = 0; j < n; j++)
		{
			b[j] = mat[j][i];
		}
		for (int j = 0; j < n; j++)
		{
			mod += b[j] * b[j];
		}
		mod = sqrt(mod);
		for (int j = 0; j < n; j++)
		{
			b[j] /= mod;
		}
		for (int j = 0; j < n; j++)
		{
			res[j][i] = b[j];
		}
	}
}
void QR(double** mat, int n, double** Q, double** R)
{
	ort_gramm_shimidt(mat, Q, n);
	norm(Q, Q, n);
	pr_mat_2(traspose(Q, n, n), mat, n, n, R);
	/*int n;
	cin >> n;
	double** mat = new double* [n];
	double** Q = new double* [n];
	double** R = new double* [n];
	for (int i = 0; i < n; i++)
	{
		mat[i] = new double[n];
		Q[i] = new double[n];
		R[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			cin >> mat[i][j];
		}
	}
	QR(mat, n, Q, R);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << Q[i][j] << "   ";
		}
		cout << endl;
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << R[i][j] << "   ";
		}
		cout << endl;
	}*/
}
void m_QR(double** mat, int n, double* lambda, double** sob_vec)
{
	double** R = new double* [n];
	double** Q = new double* [n];
	for (int i = 0; i < n; i++)
	{
		R[i] = new double[n];
		Q[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			if (i != j)
			{
				sob_vec[i][j] = 0;
			}
			else
			{
				sob_vec[i][j] = 1;
			}
		}
	}
	for (int i = 0; i < 1000; i++)
	{
		QR(mat, n, Q, R);
		pr_mat_2(traspose(Q, n, n), mat, n, n, mat);
		pr_mat_2(mat, Q, n, n, mat);
		pr_mat_2(sob_vec, Q, n, n, sob_vec);
	}
	for (int i = 0; i < n; i++)
	{
		lambda[i] = mat[i][i];
	}
}
double R(int n, int m, double(*f)(double), double a, double b)
{
	double sum = 0;
	double h;
	if (n == 0 && m == 0)
	{
		return (b - a) / 2;
	}
	else if (m == 0 && n != 0)
	{
		h = (b - a) / pow(2, n);
		for (int i = 1; i < pow(2, n - 1); i++)
		{
			sum += f(a + (2 * i - 1) * h);
		}
		return R(n - 1, 0, f, a, b) / 2 + sum * h;
	}
	else
	{
		return 1 / (pow(4, m) - 1) * (pow(4, m) * R(n, m - 1, f, a, b) - R(n - 1, m - 1, f, a, b));
	}
}
double Romberg(double(*f)(double), double a, double b)
{
	double eps = 0.000001;
	double d, res1;
	double res = R(0, 0, f, a, b);
	int i = 1;
	do {
		res1 = res;
		res = R(i, 0, f, a, b);
		d = res - res1;
		i++;
	} while (abs(d) > eps);
	return res;
}
double max_lambda(double** mat, int n)
{
	double** res = new double* [n];
	for (int i = 0; i < n; i++)
	{
		res[i] = new double[n];
	}
	double* lambda = new double[n];
	m_rot(mat, n, res, lambda);
	double max = lambda[0];
	for (int i = 0; i < n; i++)
	{
		if (max < lambda[i])
		{
			max = lambda[i];
		}
	}
	return max;
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
}
double min_lambda(double** mat, int n)
{
	double** res = new double* [n];
	for (int i = 0; i < n; i++)
	{
		res[i] = new double[n];
	}
	double* lambda = new double[n];
	m_rot(mat, n, res, lambda);
	double max = lambda[0];
	for (int i = 0; i < n; i++)
	{
		if (max > lambda[i])
		{
			max = lambda[i];
		}
	}
	return max;
}
double diff_l(double* x, double* y, int i)
{
	double h = x[i] - x[i - 1];
	double dy = y[i] - y[i - 1];
	return dy / h;
}
double diff_r(double* x, double* y, int i)
{
	double h = x[i + 1] - x[i];
	double dy = y[i + 1] - y[i];
	return dy / h;
}
double diff_c(double* x, double* y, int i)
{
	double h = x[i + 1] - x[i - 1];
	double dy = y[i + 1] - y[i - 1];
	return dy / h;
}
