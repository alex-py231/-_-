#pragma once
#include<iostream>
#include<cmath>
using namespace std;
typedef double (*f)(double* X);
double diff_l(double* x, double* y, int i);
double diff_r(double* x, double* y, int i);
double diff_c(double* x, double* y, int i);
double inter_Newton(double* X, double* Y, int n, double x);
double inter_laplas(double* X, double* Y, int n, double x);
double error_rate(double Y, double y);
double cubic_spline(double* X, double* Y, int n, double x);
void MNK(double* X, double* Y, double* res, int n, int m);
double* progonka(double** mat, int n);
double m_newtona(double (*fx)(double), double (*dfx)(double), double x0);
double m_pros_it(double (*f)(double), double a, double b);
double m_hord(double (*f)(double), double a, double b);
bool LU(double** mat, int n, int m, double** L, double** U);
double** pr_mat(double** a, double** b, int n, int m);
bool LU_clay(double** mat, double* b, double* res, int n);
//void pr_it_sis(f* f_, int n, double* res);
//void newton_sis(f* f_, f** df, int n, double* res);
void gayss(double** arr, int n, int m);
void obr(double** arr, int n);
void m_Seidel(double** mat, double* b, int n, double* res);
void m_Jacobi(double** mat, double* b, int n, double* res);
void pvr(double** mat, double* b, int n, double* res, double omega);
void m_rot(double** mat1, int n, double** res1, double* lambda);
double S_trap(double(*f)(double), double a, double b, double h);
double S_Simpson(double(*f)(double), double a, double b, double h);
double S_pr(double(*f)(double), double a, double b, double h);
void m_QR(double** mat, int n, double* lambda, double** sob_vec);
void QR(double** mat, int n, double** Q, double** R);
double Romberg(double(*f)(double), double a, double b);
double min_lambda(double** mat, int n);
double max_lambda(double** mat, int n);