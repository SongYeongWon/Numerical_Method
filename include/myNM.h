/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Song-Yeong-Won
Created          : 26-03-2018
Modified         : 29-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H
#define		PI		3.14159265358979323846264338327950288419716939937510582
#include "myMatrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//bisection and newton
double bisectionNL(double _a0, double _b0, double _tol);  // declaration of bisectionNL
double newtonRaphson(double _x0, double _tol); //declaration of newtonRaphson
double newtonRaphson_hybrid(float _x0, float _tol); // f(x)=(1/x)-2 일 때 newtonRaphson 
double hybridmethod(float _x0, float _tol); //declaration of hybridMethod
//function secant
double func1(double _x, double _y);
double funcx1d(double _x, double _y);
double funcy1d(double _x, double _y);
double func2(double _x, double _y);
double funcx2d(double _x, double _y);
double funcy2d(double _x, double _y);
// Matrix addition
extern	Matrix	addMat(Matrix _A, Matrix _B);
extern  Matrix  subMat(Matrix _A, Matrix _B); // 뺄셈 함수
extern  Matrix multMat(Matrix _A, Matrix _B); //행렬 곱 함수 
extern  void   multiMat2(Matrix _A, Matrix _B, Matrix _C);
// Apply back-substitution
extern	Matrix	backSub(Matrix _U, Matrix _b);
extern void gaussEilm(Matrix _A, Matrix _b, Matrix _U, Matrix _vecUd); // declaration of Guass elimination function
extern void LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P); //LUdecomposition 함수
extern Matrix solveLU(Matrix _L, Matrix _U, Matrix _P, Matrix _b); //solveLu 함수
extern void fwdsub(Matrix _L, Matrix _d, Matrix _y); //전진대입법 함수 
extern void inv(Matrix _A, Matrix _invA); //역행렬 구하는 함수 
extern Matrix QRdecomp(Matrix _A, Matrix _Q, Matrix _R); // QR decomposition 
extern void eigenvalue(Matrix _A, Matrix _Q, Matrix _R); // RQ function
extern double norm(Matrix _vecC);// norm 값 계산
extern Matrix transpose(Matrix _vecV); // transpose 
extern void gaussEilmJordan(Matrix _A, Matrix _b, Matrix _U, Matrix _vecUd); //guassEilmJordan

extern Matrix Fmatrix(double _x, double _y); // F matrix 반환
extern Matrix Fdmatrix(double _x, double _y); // jacobian F 프라임 matrix 반환
extern Matrix Jacobian(double _x0, double _y0); // jacobian function
extern double Cond(Matrix _A); // cond 함수


//curve fitting 
extern void linearFit(Matrix _x, Matrix _y);   // Returns the parameters of the linear least square function.
extern void planeFit(Matrix _x, Matrix _y, Matrix _z); 
extern void CurvelinearHO(Matrix _x, Matrix _y, Matrix _xq, int N);  //return interpolated values for query xq, vector or array
extern void linearInterp(Matrix _x, Matrix _y, Matrix _xq);  //return interpolated values for query xq, vector or array
extern void InterpolateLinear(Matrix _x, Matrix _y, Matrix _z);
extern void linearInterp2nd(Matrix _x, Matrix _y, Matrix _xq);
extern  Matrix arr2Mat(double* _1Darray, int _rows, int _cols); // Create a matrix from 1D-array

//Gradient Function
extern Matrix gradient(Matrix _x, Matrix _y);
extern void	gradient1D(double x[], double y[], double dydx[], int m);
extern Matrix gradientFunc(double func(const double x), Matrix xin);
extern double newtonRaphsonFunc(double func(const double x), double dfunc(const double x), double x0, double _tol);

//Integral Fucntion
extern double trapz(double x[], double y[], int m);
extern double IntegrateRect(double _x[], double _y[], int _m);
extern double integralMid(double x[], double y[], int m);
extern double integral(double func(const double x), double a, double b, int n);
extern double integral38(double func(const double x), double a, double b, int n);

//ODE Method
extern void odeEU(double func(const double x, const double y), double y[], double t0, double tf, double h,double y0);
extern void odeEM(double func(const double x, const double y), double y[], double t0, double tf, double h,double y0);
extern void odePC(double func(const double x, const double y), double y[], double t0, double tf, double h, double y0);
extern void ode(double func(const double x, const double y), double y[], double t0, double tf, double h, double	 y0,int method);
extern void odeRK2(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
extern void odeRK3(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
extern void odeRK4(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
extern void sys2RK2(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double z1[], double t0, double tf, double h, double y1_init, double z1_init);
extern void sys2RK4(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);
#endif